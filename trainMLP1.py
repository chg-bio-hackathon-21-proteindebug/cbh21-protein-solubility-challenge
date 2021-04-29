import pandas as pd
import tensorflow as tf
from tensorflow.python.keras.models import Sequential
from tensorflow.keras import layers


def train(csvFile):

	data = pd.read_csv(csvFile)

	sol_predict = data["Solubility Score"]


	data.drop(["Solubility Score","PDB File","Sequence"], inplace=True, axis=1)


	model = Sequential([
		layers.Dense(128, activation='relu', input_shape=(117,)),
		layers.Dense(128, activation='relu'),
		layers.Dense(128, activation='relu'),
		layers.Dense(64, activation='relu'),
		layers.Dense(32, activation='relu'),
		layers.Dense(16,activation='relu'),
		layers.Dense(1)
		])

	
	model.summary()


	optimizer = tf.compat.v1.keras.optimizers.Adam(learning_rate=0.001)
	
	model.compile(loss='mse',optimizer=optimizer, metrics=['mae','mse'])

	
	modelMLP = model.fit(data,sol_predict, epochs=1000, batch_size=32,verbose=1)
	print(modelMLP)
	
	tf.keras.models.save_model(model, 'model1',save_format='h5')
	
	

if __name__ == "__main__":
	
	csvFile = 'training_data.csv'
	train(csvFile)

	
