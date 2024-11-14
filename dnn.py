import os
# Keep using Keras 2
os.environ['TF_USE_LEGACY_KERAS'] = '1'
import locale
# locale.setlocale(locale.LC_ALL, 'en_US')

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

normalize = tf.keras.layers.Normalization()




loss_fn = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)

train_ds_pd, test_ds_pd = split_dataset(dataset_df)
train_features = train_ds_pd.copy()
train_labels = train_features.pop('cat')
print("{} examples in training, {} examples for testing.".format(
    len(train_ds_pd), len(test_ds_pd)))

train_features = np.array(train_features)
normalize.adapt(train_features)
model = tf.keras.Sequential([
    normalize,
    tf.keras.layers.Dense(256, activation='relu'),
    tf.keras.layers.Dense(256, activation='relu'),
    tf.keras.layers.Dense(4)
])
# Train the model.
model.compile(loss = loss_fn,
              optimizer = tf.keras.optimizers.Adam(), 
              metrics=['accuracy'])

model.fit(train_features, train_labels, epochs=20)
# evaluation = model_1.evaluate(test_ds, return_dict=True)
# print()

# for name, value in evaluation.items():
#   print(f"{name}: {value:.4f}")

# model_1.save("my_saved_model")

# tfdf.model_plotter.plot_model_in_colab(model_1, tree_idx=0, max_depth=3)