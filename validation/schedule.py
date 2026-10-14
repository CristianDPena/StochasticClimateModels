#!/usr/bin/env python3
"""Two-worker, wall-clock bounded launcher. No shell interpolation or cloud jobs."""
import argparse, csv, json, os, pathlib, shutil, signal, subprocess, time

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--batch',required=True,choices=['recovery1','recovery2','recovery3','recovery4','information1','information2','diagnostics'])
    ap.add_argument('--out',default='validation/results/laptop-v1')
    ap.add_argument('--hours',type=float,default=10.5)
    ap.add_argument('--workers',type=int,choices=[1,2],default=2)
    ap.add_argument('--max-jobs',type=int,default=0)
    args=ap.parse_args()
    if not 0 < args.hours <= 11: ap.error('--hours must be in (0,11]; leaves headroom below 12 hours')
    root=pathlib.Path(__file__).resolve().parents[1];os.chdir(root)
    out=pathlib.Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    lock=out/'scheduler.lock'
    try: lock.mkdir()
    except FileExistsError: ap.error(f'An existing scheduler lock needs inspection: {lock}')
    (lock/'pid').write_text(str(os.getpid()))
    # Calendar time includes macOS sleep; monotonic time can exclude suspended time.
    start=time.time();deadline=start+args.hours*3600;running={};log=[];prep=None
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    cmd=['Rscript','validation/run.R',f'--out={out}']
    def stop(p):
        if p.poll() is None:
            p.terminate()
            try:p.wait(timeout=20)
            except subprocess.TimeoutExpired:p.kill();p.wait()
    def clean_owned(pid):
        for d in (out/'locks').glob('*.lock'):
            if (d/'pid').exists() and (d/'pid').read_text().strip()==str(pid):shutil.rmtree(d)
    try:
        subprocess.run(cmd+['--command=init'],check=True,env=env,timeout=max(1,deadline-time.time()))
        prep_log=open(out/'logs'/f'{args.batch}_prepare.log','a')
        prep=subprocess.Popen(cmd+['--command=prepare',f'--batch={args.batch}'],stdout=prep_log,stderr=subprocess.STDOUT,env=env)
        try:code=prep.wait(timeout=max(1,deadline-time.time()))
        except subprocess.TimeoutExpired:stop(prep);return
        finally:prep_log.close()
        if code:raise RuntimeError('Cache preparation failed; inspect preparation log')
        with open(out/'design.csv',newline='') as f:rows=[r for r in csv.DictReader(f) if r['batch']==args.batch]
        rows=[r for r in rows if not (out/'results'/f"{r['job_id']}.rds").exists()]
        if args.max_jobs:rows=rows[:args.max_jobs]
        for row in rows:
            while len(running)>=args.workers and time.time()<deadline:
                for p,(r,handle,t0) in list(running.items()):
                    if p.poll() is not None:
                        handle.close();clean_owned(p.pid)
                        log.append(dict(job=r['job_id'],returncode=p.returncode,seconds=time.time()-t0))
                        del running[p]
                if len(running)>=args.workers:time.sleep(1)
            # Avoid starting a fresh costly fit with less than 45 minutes remaining.
            if time.time()+2700>=deadline:break
            handle=open(out/'logs'/f"{row['job_id']}.log",'a')
            p=subprocess.Popen(cmd+['--command=run',f"--job={row['job_id']}"],stdout=handle,stderr=subprocess.STDOUT,env=env)
            running[p]=(row,handle,time.time())
        for p,(row,handle,t0) in list(running.items()):
            try:p.wait(timeout=max(1,deadline-time.time()))
            except subprocess.TimeoutExpired:stop(p)
            handle.close();clean_owned(p.pid)
            log.append(dict(job=row['job_id'],returncode=p.returncode,seconds=time.time()-t0))
            del running[p]
        # Summary generation uses a reserved part of the <12h budget.
        subprocess.run(cmd+['--command=summarize'],check=True,env=env,timeout=900)
        subprocess.run(['Rscript','validation/figures.R',f'--out={out}'],check=True,env=env,timeout=900)
    finally:
        if prep is not None:stop(prep)
        for p,(_,handle,_) in running.items():stop(p);handle.close();clean_owned(p.pid)
        (out/'logs').mkdir(exist_ok=True)
        stamp=time.strftime('%Y%m%dT%H%M%S',time.gmtime())
        (out/'logs'/f'{args.batch}_{stamp}_scheduler.json').write_text(json.dumps(dict(batch=args.batch,elapsed=time.time()-start,jobs=log),indent=2))
        shutil.rmtree(lock)

if __name__=='__main__':main()
