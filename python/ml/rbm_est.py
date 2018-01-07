#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tensorflow as tf
import numpy as np
import os

class RBM(object):
    def __init__(self, features, num_hidden, k, persist, init_W, 
                 init_h_bias, init_v_bias, batch_size):
        self.num_vis = int(features.shape[1])
        self.num_hidden = num_hidden
        self.k = k
        self.batch_size = batch_size
        
        self.saver = None
        
        self.persistent = None
        self.persistent_update = None
        self.bit_i_idx = None
        self.bit_i_idx_update = None
        if persist:
            self.persistent = tf.Variable(tf.zeros([batch_size, num_hidden]), name="persistent")
            self.bit_i_idx = tf.Variable(0, name="bit_i_idx")
        
        self.features = features
        
        if init_W is not None:
            self.W = tf.get_variable("W", dtype=tf.float32, initializer=init_W)
        else:
            abs_val = 4*np.sqrt(6.0/(num_hidden + self.num_vis))
            self.W = tf.get_variable("W", [self.num_vis, num_hidden], dtype=tf.float32,
                                     initializer=tf.random_uniform_initializer(minval=-abs_val, maxval=abs_val))
            
        if init_h_bias is not None:
            self.h_bias = tf.get_variable("h_bias", dtype=tf.float32, initializer=init_h_bias)
        else:
            self.h_bias = tf.get_variable("h_bias", [num_hidden], dtype=tf.float32,
                                          initializer=tf.constant_initializer(0.0))
            
        if init_v_bias is not None:
            self.v_bias = tf.get_variable("v_bias", dtype=tf.float32, initializer=init_v_bias)
        else:
            self.v_bias = tf.get_variable("v_bias", [self.num_vis], dtype=tf.float32,
                                          initializer=tf.constant_initializer(0.0))
            
        self.build_model()
        
    def build_model(self):
        # Build the model:
        # We start with values for the visible nodes, perform k steps of Gibbs
        # sampling, and use the result to estimate the expected value needed by
        # the loss function

        if self.persistent is None:
            presigmoid_h, h_mean, h_sample = self.sample_h_given_v(self.features)
        else:
            h_sample = tf.stop_gradient(self.persistent)
    
        # Perform k steps of Gibbs sampling
        for _ in range(self.k):
            presigmoid_h, h_mean, h_sample, presigmoid_v, v_mean, v_sample = self.gibbs_hvh(h_sample)
        
        if self.persistent is not None:
            self.persistent_update = tf.assign(self.persistent, h_sample)
            # A hack to get the update in the loss function graph
            # TODO - Find the right way to do this
            v_sample = v_sample + tf.matmul(self.persistent_update, tf.zeros([self.num_hidden, self.num_vis]))
            
        chain_end = tf.stop_gradient(v_sample)

        self.loss = tf.reduce_mean(self.free_energy(self.features)) \
        - tf.reduce_mean(self.free_energy(chain_end))
        
        if self.persistent is None:
            self.cost = self.get_reconstruction_cost(presigmoid_v)
        else:
            self.cost = self.get_pseudo_likelihood_cost()
            
        self.saver = tf.train.Saver()
                
    def sample_h_given_v(self, v_vals):
        presigmoid_h = tf.matmul(v_vals, self.W) + self.h_bias
        h_mean = tf.sigmoid(presigmoid_h)
        h_sample = tf.nn.relu(tf.sign(h_mean - tf.random_uniform(tf.shape(h_mean))))
        return presigmoid_h, h_mean, h_sample
    
    def sample_v_given_h(self, h_vals):
        presigmoid_v = tf.matmul(h_vals,tf.transpose(self.W)) + self.v_bias
        v_mean = tf.sigmoid(presigmoid_v)
        v_sample = tf.nn.relu(tf.sign(v_mean - tf.random_uniform(tf.shape(v_mean))))
        return presigmoid_v, v_mean, v_sample
    
    def gibbs_hvh(self, h_vals):
        presigmoid_v, v_mean, v_vals = self.sample_v_given_h(h_vals)
        presigmoid_h, h_mean, h_vals = self.sample_h_given_v(v_vals)
        return presigmoid_h, h_mean, h_vals, presigmoid_v, v_mean, v_vals
    
    def gibbs_vhv(self, v_vals):
        presigmoid_h, h_mean, h_vals = self.sample_h_given_v(v_vals)
        presigmoid_v, v_mean, v_vals = self.sample_v_given_h(h_vals)
        return presigmoid_h, h_mean, h_vals, presigmoid_v, v_mean, v_vals
    
    def free_energy(self, v_vals):
        inner_prod = tf.matmul(v_vals, tf.reshape(self.v_bias, [-1,1]))
        wvpc = tf.matmul(v_vals, self.W) + self.h_bias
        sumlogs = tf.reduce_sum(tf.log(1 + tf.exp(wvpc)), reduction_indices=[1])
        return -inner_prod - tf.reshape(sumlogs, [-1,1])
    
    def get_reconstruction_cost(self, presigmoid_v):
        cross_entropy = tf.reduce_mean(tf.reduce_sum(self.features*tf.log(tf.nn.sigmoid(presigmoid_v)) + \
                                                    (1 - self.features)*tf.log(1 - tf.nn.sigmoid(presigmoid_v)), reduction_indices=[1]))
        
        #cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=self.features, logits=presigmoid_v))
        
        return cross_entropy
    
    def get_pseudo_likelihood_cost(self):
        xi = tf.round(self.features)
        fe_xi = self.free_energy(xi)
        
        i = self.bit_i_idx
        ci = tf.slice(xi, [0,i], [self.batch_size,1])
        ci_flip = 1 - ci
        xi_flip = tf.concat([tf.slice(xi,[0,0],[self.batch_size,i]), ci_flip, tf.slice(xi,[0,i+1],[self.batch_size,self.num_vis-i-1])], 1)
        fe_xi_flip = self.free_energy(xi_flip)
        
        self.bit_i_idx_update = tf.assign(self.bit_i_idx, (i + 1) % self.num_vis)
        
        return tf.reduce_mean(self.num_vis*tf.log(tf.nn.sigmoid(fe_xi_flip - fe_xi)))
    
    def train(self, sess, data_X, num_epochs, lr, beta1=0.5, checkpoint_dir=None):
        train_op = tf.train.AdamOptimizer(lr, beta1=beta1).minimize(self.loss, global_step=tf.train.get_global_step())
        sess.run(tf.global_variables_initializer())
        
        step_idx = 1
        num_batches = len(data_X) // self.batch_size
        for epoch_idx in range(num_epochs):
            loss = []
            cost = []
            for batch_idx in range(num_batches):
                #batch_xs, batch_ys = data_X.next_batch(self.batch_size)
                batch_end = (batch_idx+1)*self.batch_size
                if(batch_end > len(data_X)):
                    #batch_end = len(data_X)
                    break # TODO - Currently only acccepting exact batch sizes due to persistence
                batch_xs = data_X[batch_idx*self.batch_size:batch_end]
                
                d = {self.features:batch_xs}
                sess.run(train_op, d)
                
                loss = loss + [self.loss.eval(d, sess)]
                cost = cost + [self.cost.eval(d, sess)]
                
                # Update step
                if self.persistent is not None:
                    self.persistent_update.eval(None, sess)
                    self.bit_i_idx_update.eval(None, sess)
                
                if(np.mod(batch_idx, 100) == 0):
                    print('Batch:', batch_idx, loss[-1], cost[-1], self.bit_i_idx.eval(session=sess))
                    if checkpoint_dir is not None:
                        self.save(sess, checkpoint_dir, step_idx)
                
                step_idx = step_idx + 1
            
            print('Epoch', epoch_idx, ':', np.mean(loss), np.mean(cost))
            
            if checkpoint_dir is not None:
                self.save(sess, checkpoint_dir, step_idx)
                
    def sample(self, init_v, chain_len):
        v_sample = init_v
        for _ in range(chain_len):
            _, _, _, _, v_mean, v_sample = self.gibbs_vhv(v_sample)
            
        return v_mean, v_sample

    def save(self, sess, checkpoint_dir, step_idx):    
        model_name = "rbm.model"
    
        if not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)
        
        self.saver.save(sess, os.path.join(checkpoint_dir, model_name), global_step=step_idx)
    
    def load(self, sess, checkpoint_dir):    
        ckpt = tf.train.get_checkpoint_state(checkpoint_dir)
        if ckpt and ckpt.model_checkpoint_path:
            ckpt_name = os.path.basename(ckpt.model_checkpoint_path)
            self.saver.restore(sess, os.path.join(checkpoint_dir, ckpt_name))
            return True
        else:
            return False
    
def model_fn(features, labels, mode, params):
    
    init_W = None
    if "init_W" in params.keys():
        init_W = params["init_W"]
    
    init_h_bias = None
    if "init_h_bias" in params.keys():
        init_h_bias = params["init_h_bias"]
    
    init_v_bias = None
    if "init_v_bias" in params.keys():
        init_v_bias = params["init_v_bias"]
    
    myrbm = RBM(features["v"], params["num_hidden"], 1, False, init_W, init_h_bias, init_v_bias, params["batch_size"])
    
    train_op = tf.train.AdamOptimizer(params["learning_rate"]).minimize(myrbm.loss, global_step=tf.train.get_global_step())

    return tf.estimator.EstimatorSpec(loss=myrbm.loss, train_op=train_op, mode=mode)
    