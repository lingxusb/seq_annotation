import sklearn
import numpy as np
import matplotlib.pylab as plt
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.model_selection import cross_val_score
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.neural_network import MLPClassifier

def importdata():
    # import annotation and RNAseq reads data for both forward and backward strand
    jingm = np.loadtxt('Data\\jing_manual_cleaned.txt', dtype='int')
    M91f = np.genfromtxt(
        'D:\\Dropbox (MIT)\\Postdoc\\dataset\\Vibrio Natriegens data\\M9\\RNA-seq\\wig\\CP016345.1_f.wig')
    M91r = np.genfromtxt(
        'D:\\Dropbox (MIT)\\Postdoc\\dataset\\Vibrio Natriegens data\\M9\\RNA-seq\\wig\\CP016345.1_r.wig')
    print jingm[0, 1]

    # collect data for all promoters and terminators
    promoter = []
    for i in range(jingm.shape[0]):
        if jingm[i, 0] != 0:
            if jingm[i, 3] == 1:
                promoter.append(list(M91f[jingm[i, 0] - 50:jingm[i, 0] + 50, 1]))
            else:
                promoter.append(list(reversed(M91r[jingm[i, 0] - 50:jingm[i, 0] + 50, 1])))
    print len(promoter), promoter[0]

    term = []
    for i in range(jingm.shape[0]):
        if jingm[i, 1] != 0:
            if jingm[i, 3] == 1:
                term.append(list(M91f[jingm[i, 1] - 50:jingm[i, 1] + 50, 1]))
            else:
                term.append(list(reversed(M91r[jingm[i, 1] - 50:jingm[i, 1] + 50, 1])))
    print len(term), term[0]

    # generate random training data
    randseq = []
    while len(randseq) < 10000:
        r =  np.random.randint(60,max(jingm[:,0]))
        rexp = M91f[r-50:r+50,1]
        if np.mean(rexp) > 1:
            randseq.append(list(rexp))

    np.savetxt('Data\\jpromoters.txt', promoter,fmt='%d', delimiter='\t')
    np.savetxt('Data\\jterms.txt', term,fmt='%d', delimiter='\t')
    np.savetxt('Data\\jrandom.txt', randseq, fmt='%d', delimiter='\t')
    return 0


if __name__ == '__main__':
    #importdata()
    # load data
    jpro = np.loadtxt('Data\\jpromoters.txt', dtype='int')
    jtem = np.loadtxt('Data\\jterms.txt', dtype='int')
    jrand = np.loadtxt('Data\\jrandom.txt', dtype='int')
    print jpro[0,]

    X_train = np.concatenate((jpro[:200,], jrand[:400,]))
    Y_train = np.concatenate((np.ones((200,1)), np.zeros((400,1))))
    X = np.concatenate((jpro, jrand[:600,]))
    Y = np.concatenate((np.ones((jpro.shape[0],1)), np.zeros((600,1))))
    # method: 100 positive sample, 100 negative sample; score at CV = 5
    # SVM unoptimized results NuSVC 37/100, 0/100; SVC 39/100, 0/100, 0.82 (+/- 0.07); LinearSVC 93/100, 12/100, 0.94 (+/- 0.07)
    # svm.SVC(kernel="linear", C=0.025): 80/100, 2/100
    # Nearest neighbor: 92/100, 1/100; KNeighborsClassifier(3), 0.96 (+/- 0.04)
    # GaussianNB(): 95/100, 91/100
    # DecisionTreeClassifier(max_depth=5): 81/100, 3/100
    # RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1): 84/100, 1/100
    # AdaBoostClassifier(): 82/100, 1/100
    # Neural network: MLPClassifier(alpha=1); 91/100, 3/100, 0.85 (+/- 0.17)
    # GaussianProcessClassifier(1.0 * RBF(1.0)): 86/100,2/100
    #clf =  svm.LinearSVC(class_weight = 'balanced', C = 0.5)

    #should use LinearSVC, balanced between precision for samples and positions
    clf = svm.LinearSVC()
    scores = cross_val_score(clf, X[:,10:], Y, cv=5)
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
    clf.fit(X[:,40:], Y)
    print sum(clf.predict(X[:jpro.shape[0],10:70]))/jpro.shape[0]

    #start classification
    jingm = np.loadtxt('Data\\jing_manual_cleaned.txt', dtype='int')