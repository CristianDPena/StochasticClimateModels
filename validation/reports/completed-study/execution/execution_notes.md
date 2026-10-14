# Execution notes

The initial two-path calibration is retained in the planned study. The full
recovery1 session began on 2026-10-12 and ran across system suspensions. Long
elapsed fit times in that session are not representative CPU timings.

On 2026-10-13 at approximately 09:56 America/New_York, the first scheduler was
gracefully interrupted after 28 completed records because its monotonic deadline
could exclude macOS suspension. Its two incomplete jobs, v0149 and v0150, resumed
from existing likelihood checkpoints and subsequently completed. No completed
estimate was discarded or selectively rerun. The scheduler now uses calendar time
(commit 11d8c76); estimator and simulation source hashes are unchanged.

A sequential driver runs the seven declared batches, one at a time, with two
workers and a 10.5-hour budget per session. It stops on execution errors or a
session ending with pending jobs so these can be inspected. An idle-sleep
assertion is attached to the driver. Closing the lid or other system suspension
can still delay execution. A suspended computer enforces its deadline when it
resumes.

Elapsed times for resumed fits cover the resumed invocation and can omit work
performed before interruption. Neither these times nor suspension-affected times
should be used as unbiased method timing comparisons.
