#adapted code for YYC codec, for use in B.subtilis DNA storage project
from yyc import pipeline
from yyc import scheme
from yyc.utils import data_handle

read_file_path = "./files/teste_text.txt"
dna_path = "./output/teste_text_dna.dna"
model_path = "./output/teste_text.pkl"
write_file_path = "./output/output_teste_text.txt"

if __name__ == "__main__":
    # Primeiramente é necessário passar os parâmetros para o codec YYC. Regra Ying Yang e base 0
    [support_base, yang, ying] = ["A", #base inicial
                                    [0, 1, 0, 1], #regra Yang inicial 
                                    [[1, 1, 0, 0], [1, 0, 0, 1], [1, 1, 0, 0], [1, 1, 0, 0]]] #regra Ying inicial
    """
    Agora o codec é uma instância da classe YYC, com os parâmetros passados anteriormente
    search count: numero de iterações até aceitar o melhor scheme resultante
    max sequencia homopolimérica: 4
    max_content: porcentagem de GC
    min free energi:
    base reference: determina um padrão de bases para a codificação
    support spacing: espaçamento entre as bases de suporte para auemtar a ressitencia a erros
    """
    tool = scheme.YYC(support_bases=support_base, base_reference=yang, current_code_matrix=ying,
                      search_count=100, max_homopolymer=4, max_content=0.6)
    """
    Codificação
    Transformando o arquivo binário em sequência de DNA
    """
    pipeline.encode(
        method=tool, #modelo de codificação YCC parametrizado posteriormente
        input_path=read_file_path, 
        output_path=dna_path,
        model_path=model_path, #salva o modelo de codificação o arquivo .pkl 
        need_index=True, #indexação nas sequências para facilitar a decodificação
        segment_length=120, #tamanho do segmento de DNA
        need_log=True #mostra o log do processo
        #verify=None #método de correção de erro (aqui talvez nullomer) MAS ONDE FICA O METHOS/VERIFIES??
    )
    """
    Decodificação
    """
    del tool
    pipeline.decode(
        #method= if you want to use the same method, you can use the model file
        model_path=model_path, #ENTENDI O MODELO É SÓ OS PARAMETROS que foram passados para o YYC, podia ter a informação de index junto
        input_path=dna_path,
        output_path=write_file_path,
        has_index=True, 
        need_log=True
        # verify=None #méotod de correção de erro deveria vir aqui para fazer sentido
    )

    # compare two file
    matrix_1, _ = data_handle.read_binary_from_all(read_file_path, 120, False)
    matrix_2, _ = data_handle.read_binary_from_all(write_file_path, 120, False)
    print("source digital file == target digital file: " + str(matrix_1 == matrix_2))

    #printar o tamanho do arquivo original e do arquivo codificado em DNA
    read_file_path = open(read_file_path, "rb").read()
    print("original file size", len(read_file_path))
    dna_path = open(dna_path, "rb").read()
    print("DNA file size", len(dna_path))