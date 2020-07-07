from __future__ import print_function
from ortools.linear_solver import pywraplp
import pandas as pd
import pickle


def main():
    ## parameters
    max_peptides = 2
    affinity_threshold = 0.638
    
    
    ## load data
    path = 'C:\\Users\\Nika\\Downloads\\Organized\\Coronavirus-Epitopes\\'
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
        
    spike_covered = solver.Constraint(-1, -1, "spike_covered")
    for peptide in summary.index:
        protein = summary.at[peptide, ("Features", "protein")]
        spike_covered.SetCoefficient(x[peptide], -1 if protein in {"S1", "S2"} else 0)

    sars_cov_covered = solver.Constraint(-1, -1, "sars_cov_covered")
    for peptide in summary.index:
        covered = summary.at[peptide, ("Features", "In_SARS_Cov1")]
        sars_cov_covered.SetCoefficient(x[peptide], -1 if covered else 0)


    # objective
    objective = solver.Objective()
    for peptide in summary.index:
        df = summary.loc[peptide, "Genotypes"]
        objective.SetCoefficient(x[peptide], int(df.min()))
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
        
    # check
    for peptide in vaccine:
        print("\n")
        print(peptide)
        print("protein: " + summary.at[peptide, ("Features", "protein")])
        print("in SARS-Cov1: " + str(summary.at[peptide, ("Features", "In_SARS_Cov1")]))

if __name__ == '__main__':
    main()