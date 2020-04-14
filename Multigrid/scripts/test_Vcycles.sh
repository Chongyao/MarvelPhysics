#!/bin/bash
for gs_itrs in {1..200}
do
    ../../build/bin/V_cycle V_cycle.json $gs_itrs
done
