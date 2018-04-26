import sys
sys.path.append("../src/")
import modelSetup

if len(sys.argv)==3:
    try:
        outdir = modelSetup.outputMake(sys.argv[1],int(sys.argv[2]))
    except Exception:
        outdir = modelSetup.outputMake()
else:
    try:
        outdir = modelSetup.outputMake(sys.argv[1])
    except Exception:
        outdir = modelSetup.outputMake()
print(outdir)

