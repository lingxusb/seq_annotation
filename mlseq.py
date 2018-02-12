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
import cPickle

def import_data():
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

# harvest promoter from sequencing data
# seq: sequencing data
# chrom: DNA sequence
# model: sklearn model
# position: position of unannotated gene
# para: para[0], chrome; para[1], strand; para[2]: windwo
# pospromoter: end of the promoter region
def harvest_pro(seq, chrom, clf, position, para):
    prepromoter = []
    pospromoter = []
    seqpromoter = []
    if para[1] == -1:
        seqt = np.copy(seq[::-1])
        positiont = np.copy(position)
        positiont[:,0] = len(seqt) - position[:,1]
    else:
        seqt = np.copy(seq)
        positiont = np.copy(position)
    for i in range(positiont.shape[0]):
        if positiont[i,2] == para[0] and positiont[i,3] == para[1]:
            if seqt[positiont[i, 0],1] > 0:
                for j in range(200):
                    if clf.predict(seqt[positiont[i,0]-para[2]-j:positiont[i,0]-j,1].reshape(-1,para[2])) == 1\
                            and (seqt[positiont[i,0]-j,1]+1)/(seqt[positiont[i,0]-para[2]-j,1]+1)>2:
                        prepromoter.append(list(seqt[positiont[i,0]-para[2]-j:positiont[i,0]-j,1]))
                        seqpromoter.append(list(chrom[positiont[i,0]-100-j+para[2]:positiont[i,0]-j+para[2]]))
                        if para[1] == -1:
                            pospromoter.append([para[0], para[1], len(seqt) - positiont[i, 0] + j])
                        else:
                            pospromoter.append([para[0],para[1], positiont[i,0]-j])
                        break
    return prepromoter, pospromoter, seqpromoter

# harvest terminator from sequencing data, using the same model as promoter
# seq: sequencing data
# chrom: DNA sequence
# model: sklearn model
# position: position of unannotated gene
# para: para[0], chrome; para[1], strand; para[2]: window
# pospromoter: end of the terminator region
def harvest_ter(seq, chrom, clf, position, para):
    prepromoter = []
    pospromoter = []
    seqpromoter = []
    if para[1] == 1:
        seqt = np.copy(seq[::-1])
        positiont = np.copy(position)
        positiont[:,0] = len(seqt) - position[:,1]
    else:
        seqt = np.copy(seq)
        positiont = np.copy(position)
    for i in range(positiont.shape[0]):
        if positiont[i,2] == para[0] and positiont[i,3] == para[1]:
            if seqt[positiont[i, 0],1] > 0 and positiont[i,0]> 300:
                for j in range(200):
                    if clf.predict(seqt[positiont[i,0]-para[2]-j:positiont[i,0]-j,1].reshape(-1,para[2])) == 1\
                            and (seqt[positiont[i,0]-j,1]+1)/(seqt[positiont[i,0]-para[2]-j,1]+1)>2:
                        prepromoter.append(list(seqt[positiont[i,0]-para[2]-j:positiont[i,0]-j,1]))
                        seqpromoter.append(list(chrom[
                                                len(seqt) - positiont[i, 0] + j - 100 + para[2]:len(seqt) - positiont[
                                                    i, 0] + j + para[2]]))
                        if para[1] == 1:
                            pospromoter.append([para[0], para[1], len(seqt)-positiont[i, 0] + j])
                        else:
                            pospromoter.append([para[0],para[1], positiont[i,0]-j])
                        break
    return prepromoter, pospromoter, seqpromoter

def train_model():
    # load data
    jpro = np.loadtxt('Data\\jpromoters.txt', dtype='int')
    jtem = np.loadtxt('Data\\jterms.txt', dtype='int')
    jrand = np.loadtxt('Data\\jrandom.txt', dtype='int')
    # train models
    X_train = np.concatenate((jpro[:200,], jrand[:400,]))
    Y_train = np.concatenate((np.ones((200,1)), np.zeros((400,1))))
    neg_sample = 1000
    X = np.concatenate((jpro, jrand[:neg_sample,]))
    Y = np.concatenate((np.ones((jpro.shape[0],1)), np.zeros((neg_sample,1))))

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
    window = 40
    clf = svm.LinearSVC()
    scores = cross_val_score(clf, X[:,50-window/2:50+window/2], Y, cv=5)
    print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
    clf.fit(X[:,50-window/2:50+window/2], Y)
    #print sum(clf.predict(X[:jpro.shape[0],10:70]))/jpro.shape[0]

    # save the classifier
    with open('Data\\dumped_classifier.pkl', 'wb') as fid:
        cPickle.dump(clf, fid)
    return 0

# convert a sequence set to meme file
def output_meme():
    termseq = open('Data\\seqterminators.txt').readlines()
    f1 = open('Data\\terminators_meme.txt', 'wb')
    print termseq[0].replace(' ','')
    for i in range(len(termseq)):
        f1.write('>seq'+str(i+1))
        f1.write('\n')
        f1.write(termseq[i].replace(' ',''))
    f1.close()
    return 0

if __name__ == '__main__':
    #import_data()
    #train_model()

    # load model
    with open('Data\\dumped_classifier.pkl', 'rb') as fid:
        clf = cPickle.load(fid)

    # load sequence
    # set1: chr1 forward; set2, chr1 backward; set3 chr2 forward; set4, chr2 backward.
    # set2 and set4 are reversed
    chrom = {}
    f1 = open(
        'D:\\Dropbox (MIT)\\Postdoc\\dataset\\Vibrio Natriegens data\\index\\fasta\\Vnatriegans\\CP016345.1.dna').readlines()
    chrom["set1"] = f1[1][6:-325]
    f2 = open(
        'D:\\Dropbox (MIT)\\Postdoc\\dataset\\Vibrio Natriegens data\\index\\fasta\\Vnatriegans\\CP016346.1.dna').readlines()
    chrom["set3"] = f2[1][6:-2142]
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    chrom["set2"] = [complement.get(base, base) for base in chrom["set1"]]
    chrom["set2"] = ''.join(chrom["set2"][::-1])
    chrom["set4"] = [complement.get(base, base) for base in chrom["set3"]]
    chrom["set4"] = ''.join(chrom["set4"][::-1])

    # start classification
    unanno = np.loadtxt('Data\\unannotated.txt', dtype='int')
    M9 = {}
    M9["set1"] = np.genfromtxt(
        'D:\\Dropbox (MIT)\\Postdoc\\dataset\\Vibrio Natriegens data\\M9\\RNA-seq\\wig\\CP016345.1_f.wig')
    M9["set2"] = np.genfromtxt(
        'D:\\Dropbox (MIT)\\Postdoc\\dataset\\Vibrio Natriegens data\\M9\\RNA-seq\\wig\\CP016345.1_r.wig')
    M9["set3"] = np.genfromtxt(
        'D:\\Dropbox (MIT)\\Postdoc\\dataset\\Vibrio Natriegens data\\M9\\RNA-seq\\wig\\CP016346.1_f.wig')
    M9["set4"] = np.genfromtxt(
        'D:\\Dropbox (MIT)\\Postdoc\\dataset\\Vibrio Natriegens data\\M9\\RNA-seq\\wig\\CP016346.1_r.wig')
    # plt.plot(M9["set4"][99558-40:99558,1])
    # plt.show()
    # print chrom["set4"][-99558-60:-99558+40]

    # # classify promoters
    # window = 40
    # f1 = open('Data\\prepromoters.txt', 'ab')
    # f2 = open('Data\\pospromoters.txt', 'ab')
    # f3 = open('Data\\seqpromoters.txt', 'ab')
    # for i in range(4):
    #     prepromoter, pospromoter, seqpromoter = harvest_pro(M9["set"+str(i+1)][2:],chrom["set"+str(i+1)], clf,unanno, [i/2+1,int(((i+1)%2-0.5)*2),window])
    #     np.savetxt(f1, prepromoter, fmt='%d', delimiter='\t')
    #     np.savetxt(f2, pospromoter, fmt='%d', delimiter='\t')
    #     np.savetxt(f3, seqpromoter, fmt='%s')
    #     print len(prepromoter)
    # f1.close()
    # f2.close()
    # f3.close()
    # # classify terminators
    # f1 = open('Data\\preterminators.txt', 'ab')
    # f2 = open('Data\\posterminators.txt', 'ab')
    # f3 = open('Data\\seqterminators.txt', 'ab')
    # for i in range(4):
    #     prepromoter, pospromoter, seqpromoter = harvest_ter(M9["set"+str(i+1)][2:],chrom["set"+str(i+1)],clf,unanno, [i/2+1,int(((i+1)%2-0.5)*2),window])
    #     np.savetxt(f1, prepromoter, fmt='%d', delimiter='\t')
    #     np.savetxt(f2, pospromoter, fmt='%d', delimiter='\t')
    #     np.savetxt(f3, seqpromoter, fmt='%s')
    #     print len(prepromoter)
    # f1.close()
    # f2.close()
    # f3.close()
