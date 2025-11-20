# ROOT
This repository contains the code for the paper titled   **The origin of an irreversible transition in biological networks**

ROOT (**Revelation Of the Original circuit of irreversible Transition**) is a Python-based framework for identifying the **irreversibility kernel** underlying a specific irreversible cellular process, and for exploring **control strategies** to modulate that process.

To run ROOT, you will need:
- A Boolean network model describing the regulatory structure of the cellular system.
- Information about an irreversible phenotypic change under the irreversible cellular process.
- A mapping between external signal changes and the input configurations of the Boolean model.

Please refer to the paper for details on these inputs.

ROOT is implemented in Python 3 and is operated through a Jupyter notebook interface.

### Licence
This software is released under the MIT license.

## Installation
No special installation is required. Simply clone this repository to a suitable location as follows: 

```bash
git clone https://github.com/yena2bell/ROOT.git
```

### Requirements
- NumPy (v1.25.0) https://numpy.org/

- networkx (v3.1) https://networkx.org/

- pandas (v2.1.3) https://pandas.pydata.org/

- pyBoolnet (v3.0.11) https://github.com/hklarner/pyboolnet

## How to use ROOT
1. Open `ROOT_template/ROOT_template.ipynb` using Jupyter Notebook and ensure it's connected to a Python 3 kernel with the required packages installed.
2. In the notebook, follow the steps in order.
Cells marked with “In this cell, the user must input appropriate values” require user input.
Explanations for the required values are provided within the notebook.

### Workflow structure in the notebook:
- From the beginning of the notebook up to the section titled
"Control strategies using the irreversibility kernel",
the process focuses on identifying the irreversibility kernel.
This step is common to all ROOT applications.
- After that section, the notebook explores control strategies based on the kernel:
    - Cells before the "Resetting control" section are shared by both resetting and reversing control workflows.
    - To find resetting control strategies, execute cells from "Resetting control" up to the cell before "Reversing control".
    - To find reversing control strategies, execute cells from "Reversing control" to the end.

### Examples
You can also explore the provided examples in the `ROOT_application_examples/` folder:

These show how to apply ROOT to real biological network models.