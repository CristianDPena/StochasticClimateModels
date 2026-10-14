"""Tests the wall-clock scheduler with a fake R executable, not fake science data.
No mock output can enter the study directory.
"""
import csv, os, pathlib, subprocess, tempfile, unittest
ROOT=pathlib.Path(__file__).resolve().parents[2]
class SchedulerTests(unittest.TestCase):
    def test_budget_validation(self):
        p=subprocess.run(['python3',str(ROOT/'validation/schedule.py'),'--batch=recovery1','--hours=12'],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('11',p.stderr)
    def test_worker_limit_validation(self):
        p=subprocess.run(['python3',str(ROOT/'validation/schedule.py'),'--batch=recovery1','--workers=3'],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
    def test_tiny_budget_starts_no_fits_and_cleans_lock(self):
        with tempfile.TemporaryDirectory() as td:
            t=pathlib.Path(td);out=t/'out';exe=t/'Rscript'
            exe.write_text('''#!/usr/bin/env python3
import pathlib,sys,csv
args=dict(a[2:].split('=',1) for a in sys.argv[1:] if a.startswith('--'))
p=pathlib.Path(args['out']);(p/'logs').mkdir(parents=True,exist_ok=True);(p/'results').mkdir(exist_ok=True);(p/'locks').mkdir(exist_ok=True)
if args.get('command')=='init':
 with open(p/'design.csv','w') as f:f.write('job_id,batch\\nv0001,recovery1\\n')
if args.get('command')=='run':raise RuntimeError('Scheduler should not start an expensive fit near deadline')
''');exe.chmod(0o755)
            env=dict(os.environ,PATH=str(t)+os.pathsep+os.environ['PATH'])
            p=subprocess.run(['python3',str(ROOT/'validation/schedule.py'),'--batch=recovery1','--hours=.001',f'--out={out}'],env=env,capture_output=True,text=True,timeout=20)
            self.assertEqual(p.returncode,0,p.stderr)
            self.assertFalse((out/'scheduler.lock').exists())
            self.assertEqual(list((out/'results').glob('*.rds')),[])
if __name__=='__main__':unittest.main()
