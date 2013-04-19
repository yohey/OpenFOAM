#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import subprocess


gp = subprocess.Popen('gnuplot', shell=True,
                      stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

with gp.stdin as gin:
    gin.__write = gin.write
    gin.write = lambda s: gin.__write(bytes(s, 'utf-8'))

    gin.write('reset\n')

    gin.write('unset key\n')
    gin.write('set grid\n')
    gin.write('set logscale y\n')

    gin.write('plot "-" w lp\n')

    for n in range(5):
        with open('100/phi%d' % n, 'r') as fp:
            for line in fp:
                if 'internalField' in line:
                    for i, v in enumerate(map(float,
                                              line.strip().rstrip(');').split('(')[1].split())):
                        gin.write('%g, %g\n' % (i, v))
                    break

        gin.write('\n')
    gin.write('e\n')

    gin.flush()

    input('Press Return ...')

    gin.write('exit\n')
