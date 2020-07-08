from __future__ import print_function
from ortools.linear_solver import pywraplp
import pandas as pd
import pickle


def main():
    ## parameters
    max_peptides = 5
    affinity_threshold = 0.638

    path = 'C:\\Users\\Nika\\Downloads\\Organized\\Coronavirus-Epitopes\\'
        
    
    ## load data
    summary = pd.read_pickle(path + 'summary.pkl')
    rows = summary.index.size
    
    
    ## solver
    solver = pywraplp.Solver('simple_lp_program', pywraplp.Solver.GLOP_LINEAR_PROGRAMMING)


    # variables
    x = {}
    for peptide in summary.index:
        x[peptide] = solver.NumVar(0, 1, peptide)
        

    # constraints
    vaccine_size = solver.Constraint(0, max_peptides, "vaccine_size")
    for peptide in summary.index:
        vaccine_size.SetCoefficient(x[peptide], 1)
        
    S1_covered = solver.Constraint(-rows, -1, "S1_covered")
    for peptide in summary.index:
        protein = summary.at[peptide, ("Features", "protein")]
        S1_covered.SetCoefficient(x[peptide], -1 if protein == "S1" else 0)

    S2_covered = solver.Constraint(-rows, -1, "S2_covered")
    for peptide in summary.index:
        protein = summary.at[peptide, ("Features", "protein")]
        S2_covered.SetCoefficient(x[peptide], -1 if protein == "S2" else 0)

    E_covered = solver.Constraint(-rows, -1, "E_covered")
    for peptide in summary.index:
        protein = summary.at[peptide, ("Features", "protein")]
        E_covered.SetCoefficient(x[peptide], -1 if protein == "E" else 0)

    M_covered = solver.Constraint(-rows, -1, "M_covered")
    for peptide in summary.index:
        protein = summary.at[peptide, ("Features", "protein")]
        M_covered.SetCoefficient(x[peptide], -1 if protein == "M" else 0)

    N_covered = solver.Constraint(-rows, -1, "N_covered")
    for peptide in summary.index:
        protein = summary.at[peptide, ("Features", "protein")]
        N_covered.SetCoefficient(x[peptide], -1 if protein == "N" else 0)

    SARS_CoV1_covered = solver.Constraint(-rows, -1, "SARS_CoV1_covered")
    for peptide in summary.index:
        covered = summary.at[peptide, ("Features", "In_SARS_Cov1")]
        SARS_CoV1_covered.SetCoefficient(x[peptide], -1 if covered else 0)


    # objective
    objective = solver.Objective()
    for peptide in summary.index:
        df = summary.loc[peptide, "Genotypes"]
        objective.SetCoefficient(x[peptide], int(df.sum()))
    objective.SetMaximization()


    # solve
    solver.Solve()


    # print
    print('Number of variables =', solver.NumVariables())
    print('Number of constraints =', solver.NumConstraints())
    print('\nSolution:')
    print('Objective value =', objective.Value())
    print("\nPeptides included:")
    vaccine = set()
    for peptide in summary.index:
        included = x[peptide].solution_value() == 1
        if included:
            vaccine.add(peptide)
            print(peptide)
        
    for peptide in vaccine:
        print("\n")
        print(peptide)
        print("protein: " + summary.at[peptide, ("Features", "protein")])
        print("in SARS-Cov1: " + str(summary.at[peptide, ("Features", "In_SARS_Cov1")]))
        print("genotypes covered: " + str(summary.loc[peptide, "Genotypes"].sum()))

if __name__ == '__main__':
    main()