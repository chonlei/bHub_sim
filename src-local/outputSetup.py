import sys
sys.path.append("../src/")
import modelSetup

try:
    outdir = modelSetup.outputMake(sys.argv[1])
except Exception:
    outdir = modelSetup.outputMake()
print(outdir)

