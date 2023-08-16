# settare timout esecuzione
import os
import subprocess as sp
numberBody = 1
# 6 cicli
for i in range(1):
  print("ciao")
  print (sp.popen('ls').read())
  numberBody *= 10
  print("ciao")
  