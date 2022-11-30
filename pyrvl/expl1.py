import linalg
import op

linalg.copy('rechdr0.su','drecplh0.su')
op.fwdop('bml0.rsf','ptsrc.su','drecplh0.su')
