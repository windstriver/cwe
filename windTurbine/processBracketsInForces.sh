#!/bin/sh


cat postProcessing/turbine-forces-and-moments/0/forces.dat | tr -d '(' | tr -d ')' > postProcessing/turbine-forces-and-moments/0/forcesWithoutBrackets.dat