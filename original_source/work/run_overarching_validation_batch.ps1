param(
  [Parameter(Mandatory = $true)]
  [int]$Priority,
  [int]$MaxWorkers = 3,
  [int]$Limit = 0,
  [switch]$RetryFailed
)

$ErrorActionPreference = "Stop"
if ($MaxWorkers -lt 1 -or $MaxWorkers -gt 3) {
  throw "Frozen host limit is one to three workers."
}

$workspace = Split-Path -Parent $PSScriptRoot
$rscript = (Get-Command Rscript).Source
$runner = Join-Path $workspace "work\run_finite_window_corrected_full_pipeline.R"
$manifestPath = Join-Path $workspace "frozen_seed_manifest.csv"
$designPath = Join-Path $workspace "frozen_run_design.csv"
$hashInventoryPath = Join-Path $workspace "frozen_method_hashes.csv"
$validationRoot = Join-Path $workspace "overarching_validation_runs"
$resultsRoot = Join-Path $validationRoot "results"
$logsRoot = Join-Path $validationRoot "logs"
$recordsRoot = Join-Path $validationRoot "records"
$checkpointRoot = "D:\CodexOverarchingValidationCheckpoints"
$experimentLog = Join-Path $workspace "overarching_validation_experiment_log.csv"
$memoryLog = Join-Path $workspace "overarching_validation_memory_log.csv"

$orchestrationPaths = @(
  "work\freeze_overarching_validation_design.R",
  "work\build_overarching_correction_profiles.R",
  "work\run_overarching_validation_batch.ps1",
  "overarching_validation_protocol.md"
)
$frozenInventory = Import-Csv -LiteralPath $hashInventoryPath |
  Where-Object { $_.path -notin $orchestrationPaths }
foreach ($frozenFile in $frozenInventory) {
  $frozenPath = Join-Path $workspace $frozenFile.path
  if (-not (Test-Path -LiteralPath $frozenPath)) {
    throw "Frozen file is missing: $($frozenFile.path)"
  }
  $actualHash = (
    Get-FileHash -Algorithm SHA256 -LiteralPath $frozenPath
  ).Hash.ToLowerInvariant()
  if ($actualHash -ne $frozenFile.sha256) {
    throw "Frozen hash mismatch: $($frozenFile.path)"
  }
}

foreach ($path in @(
  $validationRoot, $resultsRoot, $logsRoot, $recordsRoot, $checkpointRoot
)) {
  New-Item -ItemType Directory -Path $path -Force | Out-Null
}

if (-not (Test-Path -LiteralPath $experimentLog)) {
  "timestamp,event,run_key,case,priority,attempt,pid,exit_code,result_path,message" |
    Set-Content -LiteralPath $experimentLog
}
if (-not (Test-Path -LiteralPath $memoryLog)) {
  "timestamp,priority,active_jobs,r_process_count,r_working_set_mb,free_physical_mb" |
    Set-Content -LiteralPath $memoryLog
}

$queue = @(
  Import-Csv -LiteralPath $manifestPath |
    Where-Object { [int]$_.priority -eq $Priority } |
    Sort-Object {[int]$_.case}
)
if ($Limit -gt 0) {
  $queue = @($queue | Select-Object -First $Limit)
}
$queue = @($queue | Where-Object {
  $result = Join-Path $resultsRoot "$($_.run_key).csv"
  $failedRecords = @(
    Get-ChildItem `
      -LiteralPath $recordsRoot `
      -Filter "$($_.run_key).attempt*.status.csv" `
      -ErrorAction SilentlyContinue |
      Where-Object {
        (Import-Csv -LiteralPath $_.FullName).status -eq "failed"
      }
  )
  (
    -not (Test-Path -LiteralPath $result) -and
    ($RetryFailed -or $failedRecords.Count -eq 0)
  )
})

$active = @{}
$nextIndex = 0
function Quote-WindowsArgument {
  param([string]$Value)
  return '"' + $Value.Replace('"', '\"') + '"'
}

function Get-ProcessTreeWorkingSetMb {
  param(
    [int]$RootProcessId,
    [object[]]$ProcessSnapshot
  )
  $treeIds = [System.Collections.Generic.HashSet[int]]::new()
  [void]$treeIds.Add($RootProcessId)
  $changed = $true
  while ($changed) {
    $changed = $false
    foreach ($item in $ProcessSnapshot) {
      if (
        $treeIds.Contains([int]$item.ParentProcessId) -and
        -not $treeIds.Contains([int]$item.ProcessId)
      ) {
        [void]$treeIds.Add([int]$item.ProcessId)
        $changed = $true
      }
    }
  }
  $workingSet = 0.0
  foreach ($item in $ProcessSnapshot) {
    if ($treeIds.Contains([int]$item.ProcessId)) {
      $workingSet += [double]$item.WorkingSetSize
    }
  }
  return $workingSet / 1MB
}

function Write-ExperimentEvent {
  param(
    [string]$Event,
    [object]$Run,
    [int]$Attempt,
    [int]$PidValue,
    [string]$ExitCodeValue,
    [string]$ResultPath,
    [string]$Message
  )
  $safeMessage = $Message.Replace('"', "'").Replace(",", ";")
  $line = @(
    (Get-Date).ToString("o"),
    $Event,
    $Run.run_key,
    $Run.case,
    $Run.priority,
    $Attempt,
    $PidValue,
    $ExitCodeValue,
    $ResultPath,
    $safeMessage
  ) -join ","
  Add-Content -LiteralPath $experimentLog -Value $line
}

while ($nextIndex -lt $queue.Count -or $active.Count -gt 0) {
  while (
    $nextIndex -lt $queue.Count -and
      $active.Count -lt $MaxWorkers
  ) {
    $run = $queue[$nextIndex]
    $nextIndex++
    $resultPath = Join-Path $resultsRoot "$($run.run_key).csv"
    $previousAttempts = @(
      Get-ChildItem `
        -LiteralPath $recordsRoot `
        -Filter "$($run.run_key).attempt*.start.csv" `
        -ErrorAction SilentlyContinue
    ).Count
    $attempt = $previousAttempts + 1
    $attemptTag = "attempt{0:D2}" -f $attempt
    $stdout = Join-Path $logsRoot "$($run.run_key).$attemptTag.out.log"
    $stderr = Join-Path $logsRoot "$($run.run_key).$attemptTag.err.log"
    $startRecord = Join-Path $recordsRoot "$($run.run_key).$attemptTag.start.csv"
    [pscustomobject]@{
      timestamp = (Get-Date).ToString("o")
      run_key = $run.run_key
      attempt = $attempt
      case = $run.case
      seed = $run.seed
      scenario_id = $run.scenario_id
      experiment = $run.experiment
      stage = $run.stage
      n_steps = $run.n_steps
      terminal = $run.terminal
      correction_profile = $run.correction_profile
      priority = $run.priority
    } | Export-Csv -LiteralPath $startRecord -NoTypeInformation

    $profilePath = Join-Path $workspace $run.correction_profile
    $argumentValues = @(
      $runner,
      "--design-file=$designPath",
      "--case=$($run.case)",
      "--tail-input=estimated",
      "--output-tag=ov_$($run.case)",
      "--chunked-likelihood",
      "--groups-per-chunk=16",
      "--kernel-cache=outputs/cache/conditional_cf_endpoint_kernel_production.rds",
      "--checkpoint-root=$checkpointRoot",
      "--correction-profile=$profilePath",
      "--validation-output=$resultPath"
    )
    if ([int64]$run.n_steps -lt 10000000) {
      $argumentValues += "--validation-resolution-experiment"
    }
    $arguments = (
      $argumentValues |
        ForEach-Object { Quote-WindowsArgument -Value $_ }
    ) -join " "
    $process = Start-Process `
      -FilePath $rscript `
      -ArgumentList $arguments `
      -WorkingDirectory $workspace `
      -WindowStyle Hidden `
      -RedirectStandardOutput $stdout `
      -RedirectStandardError $stderr `
      -PassThru
    $active[$process.Id] = [pscustomobject]@{
      Process = $process
      Run = $run
      Attempt = $attempt
      ResultPath = $resultPath
      StartTime = Get-Date
      PeakWorkingSetMb = 0
    }
    Write-ExperimentEvent `
      -Event "start" `
      -Run $run `
      -Attempt $attempt `
      -PidValue $process.Id `
      -ExitCodeValue "" `
      -ResultPath $resultPath `
      -Message "started"
  }

  Start-Sleep -Seconds 15
  $os = Get-CimInstance Win32_OperatingSystem
  $rProcesses = @(Get-Process Rscript,R -ErrorAction SilentlyContinue)
  $workingSetMb = if ($rProcesses.Count) {
    ($rProcesses | Measure-Object WorkingSet64 -Sum).Sum / 1MB
  } else {
    0
  }
  @(
    (Get-Date).ToString("o"),
    $Priority,
    $active.Count,
    $rProcesses.Count,
    [math]::Round($workingSetMb, 3),
    [math]::Round($os.FreePhysicalMemory / 1KB, 3)
  ) -join "," | Add-Content -LiteralPath $memoryLog

  $processSnapshot = @(Get-CimInstance Win32_Process)
  foreach ($pidValue in @($active.Keys)) {
    $job = $active[$pidValue]
    $job.Process.Refresh()
    if (-not $job.Process.HasExited) {
      $currentWorkingSetMb = Get-ProcessTreeWorkingSetMb `
        -RootProcessId $pidValue `
        -ProcessSnapshot $processSnapshot
      if ($currentWorkingSetMb -gt $job.PeakWorkingSetMb) {
        $job.PeakWorkingSetMb = $currentWorkingSetMb
      }
    }
    if ($job.Process.HasExited) {
      $job.Process.WaitForExit()
      $job.Process.Refresh()
      $exitCode = $job.Process.ExitCode
      $exitCodeKnown = $null -ne $exitCode
      $resultExists = Test-Path -LiteralPath $job.ResultPath
      $resultValid = $false
      if ($resultExists) {
        try {
          $resultRows = @(Import-Csv -LiteralPath $job.ResultPath)
          $resultValid = $resultRows.Count -eq 1
        } catch {
          $resultValid = $false
        }
      }
      $status = if (
        $resultValid -and
        (
          ($exitCodeKnown -and $exitCode -eq 0) -or
          -not $exitCodeKnown
        )
      ) {
        "complete"
      } else {
        "failed"
      }
      $statusPath = Join-Path $recordsRoot (
        "$($job.Run.run_key).attempt{0:D2}.status.csv" -f $job.Attempt
      )
      [pscustomobject]@{
        timestamp = (Get-Date).ToString("o")
        run_key = $job.Run.run_key
        attempt = $job.Attempt
        case = $job.Run.case
        status = $status
        exit_code = if ($exitCodeKnown) {
          $exitCode
        } else {
          "unavailable"
        }
        result_exists = $resultExists
        result_valid = $resultValid
        elapsed_seconds = ((Get-Date) - $job.StartTime).TotalSeconds
        peak_process_working_set_mb = [math]::Round(
          $job.PeakWorkingSetMb,
          3
        )
      } | Export-Csv -LiteralPath $statusPath -NoTypeInformation
      Write-ExperimentEvent `
        -Event $status `
        -Run $job.Run `
        -Attempt $job.Attempt `
        -PidValue $pidValue `
        -ExitCodeValue $(if ($exitCodeKnown) {
          [string]$exitCode
        } else {
          "unavailable"
        }) `
        -ResultPath $job.ResultPath `
        -Message $status
      $active.Remove($pidValue)
    }
  }
}

Write-Output "Priority $Priority batch complete."
