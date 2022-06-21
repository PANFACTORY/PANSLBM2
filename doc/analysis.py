import sys
import pandas as pd
import matplotlib.pyplot as plt

if(len(sys.argv) == 2):
    df = pd.read_csv(sys.argv[1])
    df_neg = df[df["htc"] < -1e-5]
    df_pos = df[df["htc"] > 1e-5]
    df_neg["htc"].hist(bins=50)
    df_pos["htc"].hist(bins=50)
    #print(df_pos["htc"].mode())
    print(df_pos["htc"].mean())
    print(df_pos["htc"].median())
    plt.xlim(-0.001, 0.005)
    plt.ylim(0, 100)
    plt.show()
    plt.close()
else:
    print("Error!")