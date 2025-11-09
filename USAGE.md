# Usage

## Environment Setup

This project uses `uv` for package management and `pyenv` to manage the Python version. The required Python version is specified in the `.python-version` file.

To set up the environment, follow these steps:

1.  **Install `pyenv` and `uv`:** Follow the official installation instructions for your operating system.
2.  **Install Python:** Navigate to the project root and run `pyenv install` to install the required Python version.
3.  **Install Dependencies:** Run `uv sync` to install the project dependencies.

## Running the Demo

To run the optimization demo, execute the following command from the project root:

```bash
/app/.venv/bin/python main.py --mode demo
```

This command runs the `main.py` script in demo mode, which will:

*   Import the necessary dependencies (`pymoo`, `pyswarms`, `matplotlib`, etc.).
*   Run the single-objective (GA and PSO) and multi-objective (NSGA-II) optimization algorithms.
*   Save the visualization plots in the `artifacts/` directory.
