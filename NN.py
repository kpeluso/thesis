import os; os.environ['TF_CPP_MIN_LOG_LEVEL']='2' # https://github.com/tensorflow/tensorflow/issues/7778
import tensorflow as tf
import numpy as np
import sys
import re

# HYPERPARAMETERS
batchSz = 20
learning_rate = 10E-4
# embedSz = 50
pi = [2.0/3.0, 1.0/3.0]

# print "Opening files..."
# trainFileName = sys.argv[1]; devFileName = sys.argv[2]
# trainFile = open(trainFileName, "r"); devFile = open(devFileName, "r");
print "Reading files..."

inpt1 = tf.placeholder(tf.int32, shape=[batchSz])
inpt2 = tf.placeholder(tf.int32, shape=[batchSz])
answr = tf.placeholder(tf.int32, shape=[batchSz])

U = tf.Variable(tf.random_normal([2*embedSz, vocabSz1], stddev=.1))
bU = tf.Variable(tf.random_normal([vocabSz1], stddev=.1))
logits = tf.matmul(both, U)+bU
# L1Output = tf.nn.relu(L1Output)
# L1_drop = tf.nn.dropout(L1Output, keep_prob)  # DROP-OUT here
# L2Output = tf.matmul(L1_drop, V)+bV
# L2Output = tf.nn.relu(L2Output)
# L2_drop = tf.nn.dropout(L2Output, keep_prob)  # DROP-OUT here

xEnt =  tf.nn.sparse_softmax_cross_entropy_with_logits(logits=logits, labels=answr)
loss = tf.reduce_mean(xEnt) # reduce mean instead of reduce sum
train = tf.train.AdamOptimizer(learning_rate).minimize(loss)

sess = tf.Session()
sess.run(tf.global_variables_initializer())
#----------------------------------------------------------
oldProg = 0
print 'Training...'
for i in range(0, tLen-batchSz, batchSz):
    next_x, next_y, next_z = trainData1[i:i+batchSz], trainData2[i:i+batchSz], trainData3[i:i+batchSz]
    l, _ = sess.run([loss, train], feed_dict={inpt1: next_x, inpt2: next_y, answr: next_z})
    progress = round(100.0*float(i)/float(tLen))
    if progress != oldProg:
    	print progress, "% done with training..."
    	oldProg = progress
print 'Testing...'
correct = 0.0; total_loss = 0.0 # reset
for i in range(0, dLen-batchSz, batchSz):
    next_x, next_y, next_z = devData1[i:i+batchSz], devData2[i:i+batchSz], devData3[i:i+batchSz]
    p = sess.run([loss], feed_dict={inpt1: next_x, inpt2: next_y, answr: next_z})
    total_loss += p[0]
    progress = round(100.0*float(i)/float(dLen))
    if progress != oldProg:
    	print progress, "% done with testing..."
    	oldProg = progress
# import pdb; pdb.set_trace()
print 'Perplexity: %.3f' % (np.exp(total_loss/(float(dLen)/float(batchSz))))

