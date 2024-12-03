import pickle

# Carrega o dicionário do arquivo pickle
with open("mona_lisa.pkl", "rb") as file:
    data_dict = pickle.load(file)

# Exibe todas as chaves e valores do dicionário
for key, value in data_dict.items():
    print("{}: {}".format(key, value))
    
    # Verifica se o valor é um objeto com atributos internos
    if hasattr(value, "__dict__"):
        print(vars(value))  # Exibe os atributos internos do objeto, se houver
