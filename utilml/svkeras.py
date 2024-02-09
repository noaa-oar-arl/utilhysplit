import datetime
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tensorflow as tf
from tensorflow import keras
import numpy as np
from keras.models import Sequential
from keras import layers
from keras.optimizers import RMSprop

# Reference - deep learning with python by Francois Chollet
# originally used to see if could predict peaks in SO2 measurements at airnow stations
# based on inputs of wind speed and direction.

def heatmap(df):
    fs = 14
    plt.matshow(df.corr(),cmap='viridis')

    plt.xticks(range(df.shape[1]), df.columns, fontsize=14,rotation=90)
    plt.gca().xaxis.tick_bottom()
    plt.yticks(range(df.shape[1]), df.columns, fontsize=14)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    plt.title("Feature Correlation Heatmap", fontsize=fs) 


def make_arma(nlen=100, ar1=[1,0.33],ma1=[1,0.9]):
    from statsmodels.tsa.arima_process import ArmaProcess
    ar1 = np.array(ar1)
    ma1 = np.array(ma1)
    sim_arma = ArmaProcess(ar1, ma1).generate_sample(nsample= nlen)
    return sim_arma

def plot_arma(sim):
    from statsmodels.tsa.stattools import pacf
    from statsmodels.tsa.stattools import acf
    from statsmodels.graphics.tsaplots import plot_pacf
    from statsmodels.graphics.tsaplots import plot_acf
    plot_pacf(sim)
    plt.show()
    plot_acf(sim)
    plt.show()


def get_hour(dt):
    return float(dt.hour)

class FakeData:

    def __init__(self):
        self.name = 'FakeData'
        self.df = pd.DataFrame()

    def create(self):
        columns1 = ['RH','TEMP','WDIR','WS','cems','hour','SO2']
        columns = ['time','RH','TEMP','WDIR','WS','cems','SO2']
        datelist = self.make_time_list()       
        nlen = len(datelist) 
        nlist = np.arange(0,nlen)
        wdir = make_arma(nlen=nlen, ar1=[1,-0.75,-0.5,0.25],ma1=[1])
        ws = make_arma(nlen=nlen, ar1=[1,-0.5,-0.2,0.10],ma1=[1])
        cems =  np.random.normal(5000,1000,nlen)
        rh = np.random.uniform(10,80,nlen)
        temp = np.random.uniform(30,50,nlen)
       
        so2 = np.zeros(nlen) 
        a1 = (wdir - wdir.mean()) / wdir.std() 
        a2 = (ws - ws.mean()) / ws.std() 
        a3 = (cems-cems.mean()) / cems.std()
        for iii in range(1,nlen):
            #so2[iii] = a1[iii-1] + a2[iii] + a3[iii]/0.5
            so2[iii] = a1[iii]*2 
        so2[0] = so2.mean()

        plist = [datelist, rh, temp,wdir,ws,cems,so2]
        plist = [np.array(x) for x in plist] 
        print(type(plist[2][0]))
        df = pd.DataFrame(plist)
        df = df.transpose()
        df.columns=columns
        df['hour'] = df.apply(lambda row: get_hour(row['time']), axis=1)
        df = df[columns1]
        for col in columns[1:]:
            df[col] = df[col].astype('float') 
        self.df = df
        return self.df

    def make_time_list(self):
        d1 = datetime.datetime(2019,1,1)
        d2 = datetime.datetime(2020,1,1)
        datelist = pd.date_range(d1, d2, freq='1H').tolist()
        return datelist 


def data_generator(data, past, future, 
                   batch_size, step, 
                   min_index=None, max_index=None, 
                   shuffle=False,target_index=-1):
    # last column of the data is the target.
    if max_index==None: 
       max_index= len(data) - future -1
    if min_index==None:
       min_index = 0
    iii=min_index + past
    done=False
    while not done:
       if iii+batch_size >=max_index:
          iii=min_index+past
       rows = np.arange(iii,min(iii+batch_size, max_index))
       iii += len(rows)
       samples = np.zeros((len(rows),
                          past // step,
                          data.shape[-1]-1))
       targets = np.zeros((len(rows),))
       for jjj, row in enumerate(rows):
          indices = range(rows[jjj]-past, rows[jjj],step)
          samples[jjj] = data[indices][:,:-1]
          targets[jjj] = data[rows[jjj] + future][target_index]
          #print(jjj, indices,rows[jjj]+future)
       yield samples, targets   


class SvLearnClass:

    def __init__(self, df,
                 batch_size=256,
                 #epochs=10,
                 split_fraction=0.5):
        self.df = df
        self.train_split = self.set_splits(split_fraction=split_fraction) 
        self.set_batch_size(batch_size)
        self.epochs=10
        #self.set_epochs(epochs)
        self.set_past(3,1)
        self.set_future(1)
        self.learning_rate = 0.001

    def set_splits(self, split_fraction=0.715):
        train_split = int(split_fraction * int(self.df.shape[0])) 
        self.train_split = train_split
        return train_split

    def set_future(self,future):
        """
        How far ahead in the future the target is.
        """
        self.future = future

    def set_past(self, past,steps):
        """
        past = number of observations in the past to be used for prediction. 
        steps = sampling rate for the observations in past
                e.g. 2 will use every other observation.
        """ 
        self.past = past
        self.steps = steps

    def set_batch_size(self, batch_size):
        self.batch_size = batch_size

    def set_epochs(self, epochs):
        self.epochs=10

    def basic(self,epoch_num=20):
        self.create_generators()
        print(self.val_steps)
        model = Sequential()
        print(self.past//self.steps, self.data.shape[-1])
        model.add(layers.Flatten(input_shape=(self.past // self.steps,
                                              self.data.shape[-1]-1)))
        model.add(layers.Dense(32, activation='relu'))
        model.add(layers.Dense(1))
        model.compile(optimizer=RMSprop(),loss='mae',metrics=['mae','acc'])
        history = model.fit(self.train_gen,
                                      steps_per_epoch=500,
                                      epochs=epoch_num,
                                      validation_data = self.val_gen,
                                      validation_steps=self.val_steps)
        self.history = history 
        self.model = model
        return model

    def gru(self,ep):
        self.create_generators()
        model = Sequential()
        print(self.past//self.steps, self.data.shape[-1])
        model.add(layers.GRU(32,input_shape=(None,
                                              self.data.shape[-1]-1)))
        #model.add(layers.Dense(32, activation='relu'))
        model.add(layers.Dense(1))
        model.compile(optimizer=RMSprop(),loss='mae')
        history = model.fit(self.train_gen,
                                      steps_per_epoch=500,
                                      epochs=ep,
                                      validation_data = self.val_gen,
                                      validation_steps=self.val_steps)
        self.history = history 
        self.model = model
        return model


    def plot_val(self):
        self.create_generators()
         


    def explore_history(self):
        sns.set()
        print(self.history.history.keys())
        loss = self.history.history['loss']
        val_loss = self.history.history['val_loss']
        epochs = range(1,len(loss)+1)
        plt.plot(epochs,loss,'bo',label='Training loss')
        plt.plot(epochs,val_loss,'b',label='Validation loss')
        plt.legend()
        plt.show()

    def predict_train(self,iii,train=True):

        if train: val = self.get_train(iii)
        else: val = self.get_val(iii)
        predicted = []
        observed = []
        for tval in zip(val[0],val[1]):
            observed.append(tval[1])
            vvv = tval[0]
            vvv = vvv.reshape(1,vvv.shape[0],vvv.shape[1])
            mval = self.model.predict(vvv)
            predicted.append(mval[0])
        return observed, predicted 

    def get_train(self,iii):
        self.create_generators()
        jjj=0
        for val in self.train_gen:
            if jjj==iii: return val
            jjj+=1
        return val

    def get_val(self,iii):
        self.create_generators()
        jjj=0
        for val in self.val_gen:
            if jjj==iii: return val
            jjj+=1
        return val

    def create_generators(self):
         
        data = self.normalize(self.df,self.train_split).values
        self.data = data
        mi = data.shape[0]-1
        self.train_gen = data_generator(data,
                                   past=self.past,
                                   future=self.future,
                                   batch_size = self.batch_size,
                                   step = self.steps,
                                   min_index=0,
                                   max_index=self.train_split)
        self.val_gen = data_generator(data,
                                   past=self.past,
                                   future=self.future,
                                   batch_size = self.batch_size,
                                   step = self.steps,
                                   min_index=self.train_split,
                                   max_index=mi)
                             
        self.val_steps = mi - self.train_split  - self.past


    def heatmap(self):
        fs = 14
        df = self.df
        plt.matshow(df.corr())

        plt.xticks(range(df.shape[1]), df.columns, fontsize=14,rotation=90)
        plt.gca().xaxis.tick_bottom()
        plt.yticks(range(df.shape[1]), df.columns, fontsize=14)
        cb = plt.colorbar()
        cb.ax.tick_params(labelsize=14)
        plt.title("Feature Correlation Heatmap", fontsize=fs) 

    def preprocess(self):
        norm_data = self.normalize(self.df, self.train_split) 
        return norm_data

    def normalize(self, df, train_split):
        data = df.values 
        data_mean = data[:train_split].mean(axis=0)
        data_std =  data[:train_split].std(axis=0)
        norm_data = (data - data_mean) / data_std
        return pd.DataFrame(norm_data)

    def train(self):
        sequence_length=24
        features = self.preprocess()
        train_data = features.loc[0:self.train_split-1]
        val_data = features.loc[train_split:]

        #x_train = train_data[[i for i range(7)]].values 
        x_train = train_data[[0,1,2,3,4,5,6]].values 
        y_train = features.iloc[start:end][[1]]

        dataset_train = keras.preprocessing.timeseries_dataset_from_array(
                              xtrain,
                              ytrain,
                              sequence_length=sequence_length,
                              sampling_rate = self.steps,
                              batch_size = self.batch_size,
                              )
 
def normalize2(df):
    data = df.values 
    data_mean = data.mean(axis=0)
    data_std =  data.std(axis=0)
    norm_data = (data - data_mean) / data_std
    return pd.DataFrame(norm_data)

def normalize(df, train_split):
    data = df.values 
    data_mean = data[:train_split].mean(axis=0)
    data_std =  data[:train_split].std(axis=0)
    norm_data = (data - data_mean) / data_std
    return pd.DataFrame(norm_data)

