import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import sys


# The code takes:
# Number of Rounds
# Sample prefix

number_of_rounds = int(sys.argv[1])
prefix = sys.argv[2]


percentage_matrix = np.empty((1000000, number_of_rounds+2), dtype=object)
i = 0
counter = 1
Rounds = []

for sample in range (0, number_of_rounds+1):
    data = pd.read_excel(io="%s_%s_Diversity.xlsx" % (prefix, sample), sheet_name="Diversity",
                           dtype={1: str, 3: float})
    Rounds = Rounds + [data]

for round in range(0, number_of_rounds+1):
    for row in range(0, len(Rounds[round])):
        print("Round is", round, "row is", row)
        sequence = Rounds[round].iloc[row][1]
        counter = round +1
        if sequence not in percentage_matrix[:,0]:
            percentage_matrix[i][0] = sequence
            percentage_matrix[i][round +1] = Rounds[round].iloc[row][3]
            for next_round in Rounds[round+1:]:
                next_round_sequence = next_round[next_round.Sequence.str.contains(sequence)]
                if len(next_round_sequence) > 0:
                    percentage_matrix[i][counter+1] = next_round_sequence.iloc[0][3]
                else:
                    percentage_matrix[i][counter+1] = 0.0
                counter +=1
            i +=1

Percentage_Columns = ["Sequence", "Original %"]
Rates_Columns = ["Sequence"]
for column in range(1, number_of_rounds+1):
    Percentage_Columns = Percentage_Columns + ["Round %s" % column]
    Rates_Columns = Rates_Columns + ["Round %s" % column]
Rates_Columns = Rates_Columns + ["Slope"]

Percentage_matrix = pd.DataFrame(percentage_matrix, columns=Percentage_Columns).dropna(thresh=2).fillna(0)

rates_matrix = np.empty((len(Percentage_matrix), number_of_rounds+2), dtype=object)
j = 0

for row in range(len(Percentage_matrix)):
    rates_matrix[row][0] = Percentage_matrix.iloc[row][0]
    for round in range(0, number_of_rounds):
        rates_matrix[row][round+1] = float(Percentage_matrix.iloc[row][round+2]) - float(Percentage_matrix.iloc[row][round+1])

x = np.arange(1, number_of_rounds+1).reshape(-1, 1)

for row in range(len(rates_matrix)):
    y = rates_matrix[row][1:number_of_rounds+1].reshape(number_of_rounds, 1)
    model = LinearRegression().fit(x, y)
    rates_matrix[row][number_of_rounds+1] = model.coef_.item()



Rates_matrix = pd.DataFrame(rates_matrix, columns=Rates_Columns).sort_values("Slope", ascending=False)

with pd.ExcelWriter("%s_Rate_analysis_%s_rounds.xlsx" % (prefix, number_of_rounds)) as writer:
    Percentage_matrix.to_excel(writer, sheet_name="Percentage", header=False)
    Rates_matrix.to_excel(writer, sheet_name="Rate", index=False)





