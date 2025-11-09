# ADME Complementarity Index (ACI) Definition

The **ADME Complementarity Index (ACI)** is a composite metric designed to replace the traditional "synergy index". It does not predict pharmacodynamic synergy, but rather evaluates the potential of a combination of traditional Chinese medicine ingredients from two dimensions: **"availability"** and **"complementarity"**, based on easily accessible physicochemical and ADME characteristics.

The optimization goal of this index is to **maximize ACI** while **minimizing hepatotoxicity**, thus forming a two-objective optimization problem.

---

## Metric Composition

ACI is composed of two weighted parts:

1.  **Mean Desirability (D)**:
    -   **Purpose**: To measure the "druggability" of each component in the combination.
    -   **Calculation**:
        1.  For each component, map its multiple ADME characteristics (such as molecular weight MW, lipophilicity AlogP, oral bioavailability OB, etc.) to an "ideal" score in the `[0, 1]` interval.
        2.  Take the **geometric mean** of these scores to obtain the final "availability" score `D_i` for that component. The geometric mean is used to ensure that an extreme shortcoming in any one characteristic will significantly lower the total score.
        3.  Finally, take the weighted (by proportion) average of the `D_i` of all components in the combination to obtain `D`.

2.  **Complementarity (C)**:
    -   **Purpose**: To measure the diversity or "complementarity" of the combination in the ADME characteristic space. A highly complementary combination means that its members have their own strengths in terms of characteristics, rather than being crowded into the same narrow area.
    -   **Calculation**:
        1.  For each ADME characteristic of all components in the combination, calculate its normalized **variance**.
        2.  Normalize the resulting variance (for example, by dividing by a preset `complement_variance_max`) to map it to the `[0, 1]` interval, to obtain the final complementarity score `C`.

---

## Calculation Formula

The calculation formula for ACI is as follows:

<p align="center">
  <strong>ACI = (1 - &lambda;) &sdot; D + &lambda; &sdot; C</strong>
</p>

-   **D**: The mean desirability of the combination.
-   **C**: The complementarity of the combination.
-   **&lambda; (lambda_weight)**: A weighting factor between `[0, 1]` used to balance the importance of availability and complementarity.
    -   When `λ` is high, combinations with large differences in characteristics among members are preferred.
    -   When `λ` is low, combinations where each member performs well on its own are preferred.
    -   The default value is usually set to `0.3`, which means that individual availability is the main factor, and complementarity is a secondary factor.

# API Documentation

This document provides detailed descriptions of the main classes and functions in `algorithms.py` and `main.py`.

---

## `algorithms.py` - Core Algorithms and Evaluation

### Data Structures

-   **`Ingredient`**: A `dataclass` used to store all relevant information for a single ingredient.
    -   **Fields**: `name`, `host_drug`, `hepatotoxicity_score`, `synergy_baseline`, and multiple ADME features (`MW`, `AlogP`, etc.).
    -   **Purpose**: Serves as the basic data unit for the entire system.

-   **`CandidateSolution`**: Represents a candidate ingredient combination.
    -   **Structure**: Internally contains a boolean list `select_bits` (determines which ingredients are selected) and a float array `proportions` (stores the proportions of each ingredient).
    -   **Main Methods**:
        -   `from_components()`: Creates an instance from specified ingredients and proportions.
        -   `iter_selects()`: Iterates over the selected ingredients and their proportions.
        -   `normalized_proportions()`: Returns the normalized proportions.

### ACI and Hepatotoxicity Evaluation Functions

These functions take a `CandidateSolution` and a list of ingredients and calculate specific performance metrics.

-   `compute_desirability_score(ingredient)`: Calculates the "desirability" score for a single ingredient.
-   `compute_complementarity_score(candidate, ingredients)`: Calculates the "complementarity" score for a combination.
-   `compute_aci_score(candidate, ...)`: Calculates the final ACI score by combining desirability and complementarity.
-   `compute_hepatotoxicity_score(candidate, ...)`: Calculates the weighted hepatotoxicity score for a combination.

### Optimization Objective Functions

These functions define the objectives for the genetic algorithm and particle swarm optimization.

-   **`single_objective_score(candidate, ...)`**: Single-objective optimization function.
    -   **Formula**: `ACI - (toxicity_weight * hepatotoxicity)`
    -   **Use**: For GA and PSO.

-   **`multi_objective_score(candidate, ...)`**: Multi-objective optimization function.
    -   **Returns**: `(-ACI, hepatotoxicity)`
    -   **Note**: Returns `-ACI` because `pymoo` minimizes all objectives by default.
    -   **Use**: For NSGA-II.

### `pymoo` and `pyswarms` Wrappers

-   **`CombinationProblem`**: A `pymoo` `Problem` subclass that connects the evaluation functions with the optimization framework.

-   **`PymooSingleObjectiveGA`**: Uses `pymoo`'s genetic algorithm to find the single-objective optimal solution.
-   **`PymooNSGAII`**: Uses `pymoo`'s NSGA-II algorithm to find the multi-objective Pareto front.
-   **`PySwarmsPSO`**: Uses `pyswarms`' particle swarm optimization algorithm to find the single-objective optimal solution.

---

## `main.py` - Program Entry Point and Flow Control

### Core Execution Functions

-   **`run_single_objective_algorithms(ingredients)`**:
    -   **Flow**:
        1.  Initializes and runs `PymooSingleObjectiveGA`.
        2.  Initializes and runs `PySwarmsPSO`.
        3.  Prints the best candidate solutions found by both algorithms.

-   **`run_nsga(ingredients, toxicity_weight)`**:
    -   **Flow**:
        1.  Initializes and runs `PymooNSGAII` to obtain multiple candidate solutions on the Pareto front.
        2.  Prints all non-dominated solutions.
        3.  Uses `select_candidate_by_weighted_sum` to select a "compromise" optimal solution from the Pareto front.

### Helper Functions

-   `sample_ingredients()`: Provides a hard-coded set of sample ingredient data for quick testing.
-   `describe_candidate(label, candidate, ingredients)`: A formatted print function to clearly display the metrics and ingredient proportions of a candidate solution.

### Program Main Entry Point

-   **`main()`**:
    1.  Calls `sample_ingredients()` to load the data.
    2.  Executes `run_single_objective_algorithms`.
    3.  Executes `run_nsga`.
