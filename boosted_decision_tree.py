import os
# Keep using Keras 2
os.environ['TF_USE_LEGACY_KERAS'] = '1'
import locale
# locale.setlocale(locale.LC_ALL, 'en_US')

import tensorflow_decision_forests as tfdf

import numpy as np
import pandas as pd
import tensorflow as tf
import tf_keras
import math

dataset_df = pd.read_csv("extracted_param_pairType_2.csv")

# Display the first 3 examples.
print(dataset_df.head(3))

label = "cat"

classes = dataset_df[label].unique().tolist()
print(f"Label classes: {classes}")

dataset_df[label] = dataset_df[label].map(classes.index)

def split_dataset(dataset, test_ratio=0.20):
  """Splits a panda dataframe in two."""
  test_indices = np.random.rand(len(dataset)) < test_ratio
  return dataset[~test_indices], dataset[test_indices]


train_ds_pd, test_ds_pd = split_dataset(dataset_df)
print("{} examples in training, {} examples for testing.".format(
    len(train_ds_pd), len(test_ds_pd)))

train_ds = tfdf.keras.pd_dataframe_to_tf_dataset(train_ds_pd, label=label)
test_ds = tfdf.keras.pd_dataframe_to_tf_dataset(test_ds_pd, label=label)


feature_1 = tfdf.keras.FeatureUsage(name="dau1_tauIdVSjet", semantic=tfdf.keras.FeatureSemantic.CATEGORICAL)
feature_2 = tfdf.keras.FeatureUsage(name="dau1_pt")
feature_3 = tfdf.keras.FeatureUsage(name="dau2_tauIdVSjet")
feature_4 = tfdf.keras.FeatureUsage(name="dau2_pt")
feature_5 = tfdf.keras.FeatureUsage(name="pNet_sum")
feature_6 = tfdf.keras.FeatureUsage(name="bjet1_pt")
feature_7 = tfdf.keras.FeatureUsage(name="bjet2_pt")
feature_8 = tfdf.keras.FeatureUsage(name="bjet1_eta")
feature_9 = tfdf.keras.FeatureUsage(name="bjet2_eta")
all_features = [feature_1, feature_2, feature_3, feature_4, feature_5, feature_6, feature_7, feature_8, feature_9]

# Specify the model.
model_1 = tfdf.keras.GradientBoostedTreesModel(growing_strategy="LOCAL", num_trees=1000, max_depth=5,
    split_axis="SPARSE_OBLIQUE",
    categorical_algorithm="RANDOM",)

# Train the model.
model_1.fit(train_ds)

model_1.compile(metrics=["accuracy"])
evaluation = model_1.evaluate(test_ds, return_dict=True)
print()

for name, value in evaluation.items():
  print(f"{name}: {value:.4f}")

model_1.save("my_saved_model")

tfdf.model_plotter.plot_model_in_colab(model_1, tree_idx=0, max_depth=3)