#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tensorflow as tf
from tensorflow.python import debug as tf_debug
import numpy as np
import rbm_est
from tensorflow.python import debug as tf_debug

from tensorflow.examples.tutorials.mnist import input_data
mnist = input_data.read_data_sets("/Users/tuckerkj/data/MNIST_data/", one_hot=True)

tf.set_random_seed(1)

batch_size = 64
k = 15
checkpoint_dir = '/Users/tuckerkj/checkpoints/rbm_mnist_pers2_k15_b64'

train_input_fn = tf.estimator.inputs.numpy_input_fn(x={"v":mnist.train.images},
                                                    y=mnist.train.labels,
                                                    batch_size=batch_size,
                                                    num_epochs=25,
                                                    shuffle=True)

params = {"k":k, "learning_rate":0.01, "num_hidden":500, "num_vis":784, "batch_size":batch_size}
#rbm = tf.estimator.Estimator(model_fn=rbm_est.model_fn, params=params, model_dir='/tmp/rbm3')

tf.logging.set_verbosity(tf.logging.INFO)

#rbm.train(input_fn=train_input_fn, steps=2000)

features = tf.placeholder(tf.float32, [None, 784], "features")
myrbm = rbm_est.RBM(features, 500, k, True, None, None, None, batch_size)

sess = tf.Session()
#with tf.Session() as sess:
    #sess = tf_debug.LocalCLIDebugWrapperSession(sess)
myrbm.train(sess, mnist.train.images, 15, 0.01, checkpoint_dir=checkpoint_dir)