import os
import joblib
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from moment_matching_evaluation import *

def plot_grid_scores(data):
    a00vals = data["a00vals"]
    a10vals = data["a10vals"]
    results = data["results"]
    
    score_grid = np.zeros((len(a10vals), len(a00vals)))
    
    idx = 0
    for i in range(len(a10vals)):
        for j in range(len(a00vals)):
            res = results[idx]
            # Average the random forest scores across the n_average attempts
            scores = [attempt["score"] for attempt in res["tries"]]
            score_grid[i, j] = np.mean(scores)
            idx += 1

    plt.figure(figsize=(8, 6))
    # Round ticks for clean display
    xticklabels = [f"{x:.1f}" for x in a00vals]
    yticklabels = [f"{y:.1f}" for y in a10vals]
    
    sns.heatmap(
        score_grid, 
        xticklabels=xticklabels, 
        yticklabels=yticklabels, 
        annot=True, 
        fmt=".2f", 
        cmap="viridis",
        cbar_kws={'label': 'Random Forest Accuracy'}
    )
    plt.xlabel(r"Variance of $x_0$ ($a_{00}$)")
    plt.ylabel(r"Coupling Norm ($a_{i0}$)")
    plt.title("Classifier Separation Score Across Parameter Grid")
    plt.gca().invert_yaxis() # Keep lower values at the bottom-left
    plt.tight_layout()

def plot_samples_at_coordinate(data, a00_target, ai0_target, attempt_idx=0):
    """
    Finds the closest grid coordinates and plots the underlying distribution samples 
    comparing Wrapped Normal vs Generalized Von Mises.
    """
    results = data["results"]
    
    best_match = None
    min_dist = float('inf')
    
    for res in results:
        dist = (res["a00"] - a00_target)**2 + (res["ai0"] - ai0_target)**2
        if dist < min_dist:
            min_dist = dist
            best_match = res
            
    if best_match is None:
        print("No matching grid coordinate found.")
        return

    # Grab the requested attempt index data
    attempt_data = best_match["tries"][attempt_idx]
    norm_samples = attempt_data["norm_samples"]
    gvm_samples = attempt_data["gvm_samples"]
    
    num_samples, total_dim = norm_samples.shape
    
    print(f"Plotting samples for closest grid point: a00={best_match['a00']:.2f}, ai0={best_match['ai0']:.2f}")
    print(f"Total Dimensions detected: {total_dim} (1 Angular + {total_dim - 1} Euclidean)")

    # Prepare data for Seaborn (melt into long format or stack)
    labels = np.array(["Wrapped Normal"] * num_samples + ["GVM"] * num_samples)
    combined_data = np.vstack((norm_samples, gvm_samples))
    
    # Define Column Names dynamically
    col_names = ["Angular ($x_0$)"] + [f"Euclidean ($x_{i}$)" for i in range(1, total_dim)]
    
    # Create master visualization path based on dimensions
    if total_dim == 2:
        # Simple 2D Scatter plot
        plt.figure(figsize=(7, 5))
        sns.scatterplot(
            x=combined_data[:, 0], 
            y=combined_data[:, 1], 
            hue=labels, 
            alpha=0.1, 
            palette="Set1",
            edgecolor=None
        )
        plt.xlabel(col_names[0])
        plt.ylabel(col_names[1])
        plt.title(f"Distribution Comparison\n($a_{{00}}$={best_match['a00']:.1f}, $a_{{i0}}$={best_match['ai0']:.1f})")
        plt.tight_layout()
        
    elif total_dim == 3:
        # 3D scatter plot via matplotlib (seaborn doesn't natively do 3D projections)
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        
        for lbl, color in [("Wrapped Normal", "red"), ("GVM", "blue")]:
            mask = (labels == lbl)
            ax.scatter(
                combined_data[mask, 0], 
                combined_data[mask, 1], 
                combined_data[mask, 2], 
                alpha=0.4, 
                label=lbl,
                edgecolor=None,
                s=15
            )
        ax.set_xlabel(col_names[0])
        ax.set_ylabel(col_names[1])
        ax.set_zlabel(col_names[2])
        plt.title(f"3D Distribution Comparison\n($a_{{00}}$={best_match['a00']:.1f}, $a_{{i0}}$={best_match['ai0']:.1f})")
        plt.legend()
        plt.tight_layout()
        
    else:
        # High Dimension Handler: Vectorized Pair Plot Matrix
        import pandas as pd
        print("Dimension > 3. Generating a Seaborn Pair Plot grid...")
        df = pd.DataFrame(combined_data, columns=col_names)
        df["Distribution"] = labels
        
        # Build the pairwise correlation grid plot
        g = sns.pairplot(
            df, 
            hue="Distribution", 
            palette="Set1", 
            plot_kws={'alpha': 0.3, 'edgecolor': None, 's': 10},
            diag_kind="kde"
        )
        g.fig.suptitle(f"Pairwise Feature Projections ($Dim={total_dim}$)", y=1.02)


if __name__ == "__main__":
    sweep_data = joblib.load("results.pkl")
    
    plot_grid_scores(sweep_data)
    
    target_a00 = 1.1
    target_ai0 = 1.5
    plot_samples_at_coordinate(sweep_data, a00_target=target_a00, ai0_target=target_ai0, attempt_idx=0)
    
    plt.show()
