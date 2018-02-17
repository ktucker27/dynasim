#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tensorflow as tf
from tensorflow.python import debug as tf_debug
import numpy as np
import rbm_est
from tensorflow.python import debug as tf_debug

from tensorflow.examples.tutorials.mnist import input_data
mnist = input_data.read_data_sets("/Users/tuckerkj/data/MNIST_data/", one_hot=True)

NUM_THREADS = 4

tf.set_random_seed(1)

batch_size = 64
k = 15
checkpoint_dir = '/Users/tuckerkj/checkpoints/rbm_mnist_pers2_k15_b64'
n_samples = 10
n_chains = 20
chain_len = 1000

features = tf.placeholder(tf.float32, [None, 784], "features")
myrbm = rbm_est.RBM(features, 500, k, True, None, None, None, batch_size)

sess = tf.Session(config=tf.ConfigProto(intra_op_parallelism_threads=NUM_THREADS))
myrbm.load(sess, checkpoint_dir)

start_idx = 100
init_v_vals = mnist.test.images[start_idx:start_idx+n_chains]
np.savetxt("/Users/tuckerkj/data/scratch/sample_in.csv", init_v_vals, delimiter=",")
init_v = tf.placeholder(tf.float32, [None, 784], "init_v")
all_v_means = None
for i in range(n_samples):
    print('Sample:', i)
    v_mean, v_sample = myrbm.sample(init_v, chain_len)
    v_mean_vals, init_v_vals = sess.run([v_mean, v_sample], feed_dict={init_v:init_v_vals})
    if all_v_means is None:
        all_v_means = v_mean_vals
    else:
        all_v_means = np.append(all_v_means, v_mean_vals, axis=0)

    np.savetxt("/Users/tuckerkj/data/scratch/sample_out2.csv", all_v_means, delimiter=",")