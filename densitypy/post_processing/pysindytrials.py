import gc

import numpy as np
import pandas as pd
import pysindy as ps
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy

from densitypy.project_utils.daskloading import fetch_csv_data_dask


def dummy_contraints():
    from pysindy.optimizers import ConstrainedSR3
    from pysindy.feature_library import PolynomialLibrary
    from pysindy import SINDy

    # Define your constraints matrix A and vector b (A*X = b)
    # This is a simple example. Replace it with your actual constraints.
    constraint_lhs = [[0, 1, -1, 0], [0, 0, 1,
                                      -1]]  # Example constraint: second and third features are equal, third and fourth are equal
    constraint_rhs = [0, 0]  # The equations equal zero

    # Initialize the ConstrainedSR3 optimizer with your constraints
    optimizer = ConstrainedSR3(constraint_lhs=constraint_lhs, constraint_rhs=constraint_rhs, threshold=0.1)

    # Initialize SINDy as before, but with the new optimizer
    model = SINDy(optimizer=optimizer, feature_library=PolynomialLibrary(degree=2))

    # Fit and print the model as before
    model.fit(X, t=t)
    model.print()


def fetch_atomic_charge_all_csv(experiment_path, columns_to_drop=None, isim=189):
    TIME_COLUMN_NAME = 'itime'

    atomic_charege_file_path = f'{experiment_path}/AtomicCharge/AtomicCharge_ALL.csv'

    # ------------------------LOAD DATA-----------------------------
    observable = fetch_csv_data_dask(data_file_path=atomic_charege_file_path)

    # ---------------------USER CUSTOM DATA PREPROCESSING----------------
    if columns_to_drop is not None:
        observable = observable.drop(columns=columns_to_drop)

    observable = observable[observable['iSim'] == isim]

    # ---------------------Pysindy Required DATA PREPROCESSING----------------
    observable.sort_values(by=TIME_COLUMN_NAME, inplace=True)

    t = observable.pop(TIME_COLUMN_NAME).values
    X = observable

    gc.collect()
    return X, t


def prep_for_atomic_charge_sim_csv(experiment_path,
                                   column_names_to_drop=None,
                                   time_column_name: str = "Time"):
    """

    # "itime","Time","TotalCharge","Atom_O_ChargeX","Atom_O_ChargeY","Atom_O_ChargeZ",
    # "Atom_N_ChargeX","Atom_N_ChargeY","Atom_N_ChargeZ", ...

    :param experiment_path:
    :param column_names_to_drop:
    :param time_column_name:
    :return:
    """
    if column_names_to_drop is None:  # This is a hack to avoid mutable default arguments
        column_names_to_drop = ["itime", 'TotalCharge']

    atomic_charege_file_path = f'{experiment_path}/AtomicCharge/AtomicChargeXUV.csv'

    # ------------------------LOAD DATA-----------------------------
    observable = fetch_csv_data_dask(data_file_path=atomic_charege_file_path)

    # ---------------------USER CUSTOM DATA PREPROCESSING (Will be scratches for now)----------------
    if column_names_to_drop is not None:
        observable = observable.drop(columns=column_names_to_drop)

    t = observable.pop(time_column_name).values
    X = observable

    gc.collect()
    return X, t


def run_pipeline(X, t):
    # Initialize a differentiation method
    # differentiation_method = ps.differentiation.FiniteDifference(order=3)
    differentiation_method = ps.differentiation.SpectralDerivative(d=3)

    # Initialize the feature library.
    library_features_list = [
        ps.feature_library.FourierLibrary(
            n_frequencies=4,
            include_sin=True,
            include_cos=True,
            library_ensemble=True),
        ps.feature_library.PDELibrary(),
        # ps.PolynomialLibrary(degree=2),
    ]
    library_features = None
    for lib in library_features_list:
        if library_features is None:
            library_features = lib
        else:
            library_features = library_features * lib
    library_features_list.append(library_features)

    feature_library = ps.feature_library.GeneralizedLibrary(library_features_list)
    feature_library.fit(X)
    feature_library.transform(X)

    # feature_library = None
    # for library_features in library_features_list:
    #     if feature_library is not None:
    #         feature_library = feature_library + library_features
    #     else:
    #         feature_library = library_features

    # Initialize the optimizer, e.g., STLSQ optimizer
    optimizer = ps.STLSQ(
        threshold=0.01,
        alpha=0.05,
        max_iter=20,
        ridge_kw=None,
        normalize_columns=True,
        fit_intercept=False,
        copy_X=True,
        initial_guess=None,
        verbose=True
    )
    optimizer = ps.optimizers.SINDyOptimizer(optimizer=optimizer, unbias=False)
    # optimizer = ps.EnsembleOptimizer(opt=optimizer, library_ensemble=True)

    # Initialize and fit the SINDy model focusing on the 'close' column dynamics
    model = ps.SINDy(
        differentiation_method=differentiation_method,
        feature_library=feature_library,
        optimizer=optimizer,
        feature_names=X.columns,
    )
    model.fit(X, t=t)

    model.print()

    # Predict derivatives using the learned model
    x_dot_test_predicted = model.predict(X)

    # Compute derivatives with a finite difference method, for comparison
    x_dot_test_computed = model.differentiate(X)

    # Dont plot just compare predicted with computed and provide a mean_squared_error
    RMSE = ((x_dot_test_predicted - x_dot_test_computed) ** 2).mean() ** 0.5
    print(f'Root Mean Squared Error: {RMSE}')

    # fig, axes = plt.subplots(nrows=X.shape[1], ncols=1, figsize=(12, 12))
    # for i in range(X.shape[1]):
    #     print(f'Plotting {i} of {X.shape[1]}')
    #     axes[i].plot(t, x_dot_test_computed[:, i], "k", label="numerical derivative")
    #     axes[i].plot(t, x_dot_test_predicted[:, i], "r--", label="model prediction")
    #     axes[i].set(xlabel="t", ylabel=r"$\dot x_{}$".format(i))
    # plt.tight_layout()
    # plt.show()


def optimal_hierarchical_cluster(mat: pd.DataFrame, method: str = "ward") -> np.array:
    """
    Calculates the optimal clustering of a matrix.

    It calculates the hierarchy clusters from the distance of the matrix. Then it calculates
    the optimal leaf ordering of the hierarchy clusters, and returns the optimally clustered matrix.

    It is reproduced with modifications from the following blog post:
    `Marti, G. (2020) TF 2.0 DCGAN for 100x100 financial correlation matrices [Online].
    Available at: https://marti.ai/ml/2019/10/13/tf-dcgan-financial-correlation-matrices.html.
    (Accessed: 17 Aug 2020)
    <https://marti.ai/ml/2019/10/13/tf-dcgan-financial-correlation-matrices.html>`_

    This method relies and acts as a wrapper for the `scipy.cluster.hierarchy` module.
    `<https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html>`_

    :param mat: (np.array/pd.DataFrame) Correlation matrix.
    :param method: (str) Method to calculate the hierarchy clusters. Can take the values
        ["single", "complete", "average", "weighted", "centroid", "median", "ward"].
    :return: (np.array) Optimal hierarchy cluster matrix.
    """

    matrix = np.array(mat)
    labels = np.array(mat.columns)
    d = hierarchy.distance.pdist(matrix)  # vector of upper triangle of distance matrix
    z = hierarchy.linkage(d, method=method)  # linkage matrix
    optimal_leaf_ordering = hierarchy.leaves_list(hierarchy.optimal_leaf_ordering(z, d))
    optimal_hierarchy_cluster = pd.DataFrame(np.array(matrix[optimal_leaf_ordering, :])[:, optimal_leaf_ordering],
                                             index=labels[optimal_leaf_ordering],
                                             columns=labels[optimal_leaf_ordering])

    return optimal_hierarchy_cluster


def plot_optimal_hierarchical_cluster(mat, method="ward"):
    """
    Calculates and plots the optimal clustering of a correlation matrix.

    It uses the `optimal_hierarchical_cluster` function in the clustering module to calculate
    the optimal hierarchy cluster matrix.

    :param mat: (np.array/pd.DataFrame) Correlation matrix.
    :param method: (str) Method to calculate the hierarchy clusters. Can take the values
        ["single", "complete", "average", "weighted", "centroid", "median", "ward"].
    :return: (plt.Axes) Figure's axes.
    """

    # Calculate optimal hierarchical cluster
    cluster = optimal_hierarchical_cluster(mat, method)

    # Plot optimal hierarchical cluster
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))
    ax = sns.heatmap(cluster, cmap="RdBu", annot=True, fmt=".2f", ax=ax, center=0)
    ax.set_title("Optimal Hierarchical Cluster")
    ax.set_xlabel("Asset")
    ax.set_ylabel("Asset")

    return ax


if __name__ == "__main__":
    EXPERIMENT_PATH = "/home/ruben/PycharmProjects/DensityPy/Studies/ExampleStudy/debug_sim"
    # X, t = prep_for_atomic_charge_sim_csv(EXPERIMENT_PATH)
    # correlation = X.corr().dropna(axis=0, how='all').dropna(axis=1, how='all')
    # columns_not_in_correlation = [column_name for column_name in X.columns if column_name not in correlation.columns]
    #
    # for column_name in columns_not_in_correlation:
    #     X.drop(column_name, axis=1)
    #
    # plt.figure(figsize=(6, 4))
    # plt.pcolormesh(correlation, cmap='RdBu', vmin=-1, vmax=1)
    # # Ticks and labels
    # plt.xticks(np.arange(0, len(correlation.columns), 1), correlation.columns, rotation=90)
    # plt.yticks(np.arange(0, len(correlation.index), 1), correlation.index)
    # plt.colorbar()
    # plt.title("Original Correlation Matrix")
    #
    # # Plot optimal clustering.
    # plot_optimal_hierarchical_cluster(correlation, method="ward")
    # plt.title("Optimal Clustering Correlation Matrix")
    # plt.show()

    X, t = fetch_atomic_charge_all_csv(EXPERIMENT_PATH)
    run_pipeline(X, t)
