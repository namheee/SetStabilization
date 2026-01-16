# SetStabilization

## Usage

+ You can simply download SetStabilization from this git repository, while setup.py is not provided.
+ We provide a main simulation example using Jupyter notebooks. Please refer to `Tutorial.ipynb`.
+ An example network is provided in the `example` directory.  


## Notes
+ To generate random Boolean networks, we utilized [BNGenerator](https://github.com/choonlog/OutputStabilization), a software tool that constructs random networks based on biological Boolean logics extracted from the Cell Collective (https://cellcollective.org/).
+ To identify all attractors of Boolean networks, we utilized [BooleanSim](https://github.com/jehoons/BooleanSim), a Python 3 tool developed based upon the [booleannet](https://github.com/ialbert/booleannet).
+ Main functions are included in `setMain.py` and `setFunction.py`.
+ An auxiliary function is included in `choonFunction.py`, with a referecne to its source at [AttractorTransition](https://github.com/choonlog/AttractorTransition).


## Comparison with other algorithms
+ For brute force search and LDOI control, we utilited [pystablemotifs](https://github.com/jcrozum/pystablemotifs).
```python
pystablemotifs.drivers.minimal_drivers for brute force search

pystablemotifs.drivers.GRASP for greedy randomized adaptive search procedure (GRASP) search (GRASP_iterations = 2,000) for LDOI control
```

+ For [stable motif control](https://github.com/jcrozum/pystablemotifs),
```python
ar=pystablemotifs.AttractorRepertorie.from_primes and ar.succession_diagram.reprogram_to_trap_spaces(); specifically, ‘AttractorRepertoire.from_primes’ is employed to identify the complete set of attractors, while ‘reprogram_to_trap_spaces’ is used to determine the control inputs based on the succession diagram; the computational time includes all steps, including attractor identification and control input search.
```

+ For [FVS control](https://github.com/CASCI-lab/CANA),
```python
cana.control.feedback_vertex_set_driver_nodes(graph='structural', method='bruteforce', max_search=10, keep_self_loops=True)
```

+ For [biobalm control](https://github.com/jcrozum/biobalm)
```python
biobalm.SuccessionDiagram.from_rules and biobalm.control.succession_control()
```


## Reference paper
+ A related paper is expected to be published soon.
