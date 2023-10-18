#!/usr/bin/env python

'''
plotting_functions_misc_akp.py

Helper functions for plotting performance of models in HP searches

DILI transporter classification and regression models
'''


import sklearn
from sklearn.metrics import confusion_matrix
import itertools
import matplotlib.pyplot as plt
import numpy as np

# confusion matrix
# from sklearn 0.19.2 documentation:
def plot_confusion_matrix(cm, classes,
                          normalize=False,
                          cmap=plt.cm.Blues, ax=None):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
#         print("Normalized confusion matrix")
#     else:
#         print('Confusion matrix, without normalization')
#     print(cm)
    if ax is None:
        im = plt.imshow(cm, interpolation='nearest', cmap=cmap)
#         plt.colorbar(cmap=cmap, shrink=0.7)
        tick_marks = range(0,len(classes))
        plt.xticks(tick_marks, classes, rotation=0)
        plt.yticks(tick_marks, classes, rotation=90, va='center')

        fmt = '.2f' if normalize else 'd'
        thresh = cm.max() / 2.
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            plt.text(j, i, format(cm[i, j], fmt),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")
        plt.tight_layout()
        plt.ylabel('True label')
        plt.xlabel('Predicted label')
    else:
        im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
#         plt.colorbar(im, shrink=0.7)
        fmt = '.2f' if normalize else 'd'
        thresh = cm.max() / 2.
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            ax.text(j, i, format(cm[i, j], fmt),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")
#         ax.tight_layout()
        ax.set_ylabel('True label')
        ax.set_xlabel('Predicted label')
        tick_marks = range(0,len(classes))
        ax.set_xticks(tick_marks)
        ax.set_xticklabels(classes, rotation=0)
        ax.set_yticks(tick_marks)
        ax.set_yticklabels(classes, rotation=90, va='center')
    return im