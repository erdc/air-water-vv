Rubble-Mound Breakwater
==============================================

Running the test case
-----

To run in parallel (example with mpirun and 8 processors):

```
mpirun -np 8 parun --TwoPhaseFlow crownwall_breakwater.py -D test2 -v -C "mwl=0.671 Hs=0.158 Tp=1.855"
```

```
mpirun -np 8 parun --TwoPhaseFlow crownwall_breakwater.py -D test1 -v -C "mwl=0.671 Hs=0.190 Tp=2.032"
```


```
mpirun -np 8 parun --TwoPhaseFlow crownwall_breakwater.py -D test3 -v -C "mwl=0.671 Hs=0.145 Tp=1.721"
```

