"""
Esse código é um script python para um projeto para armazenamento de dados em DNA em B.subtilis
Arthur Busanello UFRGS 2024
"""


# A coisa mais importante
import sys
sys.path.insert(0, "/home/busa/projetos/dna-storage/yyc/DNA-storage-YYC-available")
## Importando as bibliotecas necessárias
### Ying Yang Codec
from yyc import pipeline
from yyc import scheme
### Outras bibliotecas
import os
import pandas as pd

## Inicialização de paths de arquivos e diretórios do projeto
input_file_path = "busa-inputs/textual.txt"
input_partes_dir = "busa-inputs/partes"
os.makedirs(input_partes_dir, exist_ok=True) #aqui vai ser as partes do texto tratado
output_dir = "busa-outputs"
output_partes_dir = os.path.join(output_dir, "partes") #aqui vai ser as partes codificadas
os.makedirs(output_partes_dir, exist_ok=True)
dna_dir = os.path.join(output_dir, "dna")
os.makedirs(dna_dir, exist_ok=True)
decoded_dir = os.path.join(output_dir, "decoded")
os.makedirs(decoded_dir, exist_ok=True)
genome_path = "busa-inputs/Bsub-Cohn-genome.fasta"

## Parâmetros para o tratamento do arquivo de entrada
tamanho_parte = 1200 #tamanho que cada parte do arquivo de entrada vai ser dividido
segment_length = 120 #tamanho do segmento dos segmentos de DNA da codificação



## Main 
if __name__ == "__main__":

    ## Tratamento inicial do arquivo de entrada
    print("Tratamento inicial do arquivo de entrada...")
    with open(input_file_path, "r") as file:
        data = file.read()
    
    partes = [data[i:i+tamanho_parte] for i in range(0, len(data), tamanho_parte)]
    
    for idx, parte in enumerate(partes):
        parte_path = os.path.join( input_partes_dir, "parte_{}.txt".format(idx+1))
        with open(parte_path, "w") as parte_file:
            parte_file.write(parte)

    ### log partes divididas
    for idx, parte in enumerate(partes):
        print("Parte {}: Tamanho = {} caracteres".format(idx+1, len(parte)))


    ## Parâmetros do Ying Yang Codec
    ### Inicialização do codec, base de suporte, regras Ying e Yang
    [support_base, yang, ying] = ["A", #base inicial
                                    [0, 1, 0, 1], #regra Yang inicial 
                                    [[1, 1, 0, 0], [1, 0, 0, 1], [1, 1, 0, 0], [1, 1, 0, 0]]] #regra Ying inicial
    
    ### Instância do codec YYC
    tool = scheme.YYC(support_bases=support_base, 
                      base_reference=yang, 
                      current_code_matrix=ying,
                      search_count=100, 
                      max_homopolymer=4, 
                      max_content=0.6
                      # min_free_energy= None
                      )
    


        ## Codificação / Enconding
    print("CODIFICAÇÃO DAS PARTES DO ARQUIVO DE ENTRADA EM DNA...")
    ### Transformando o arquivo binário em sequência de DNA
    data_final = []
    primer_5 = "ATCG"  # Primer 5' para a sequência de DNA
    primer_3 = "GCTA"  # Primer 3' para a sequência de DNA
    # Caminho do log consolidado
    log_path = os.path.join(output_dir, "log_codificacao.md")
    with open(log_path, "w") as log_file:
        log_file.write("# Log de Codificação e Decodificação\n\n")

    # Lista para armazenar os dados da tabela
    tabela_dados = []

    for idx, parte_path in enumerate(sorted(os.listdir(input_partes_dir), key=lambda x: int(x.split("_")[-1].split(".")[0]))):
        input_path = os.path.join(input_partes_dir, parte_path)
        dna_output_path = os.path.join(output_partes_dir, "parte_{}.dna".format(idx + 1))
        model_path = os.path.join(output_partes_dir, "parte_{}.pkl".format(idx + 1))

        print("Codificando Parte {}: {} -> {}".format(idx + 1, input_path, dna_output_path))

        pipeline.encode(
            method=tool,  # modelo de codificação YCC parametrizado
            input_path=input_path,
            output_path=dna_output_path,  # salva a parte codificada em DNA
            model_path=model_path,  # salva o modelo de codificação o arquivo .pkl
            need_index=True,  # indexação nas sequências para facilitar a decodificação
            segment_length=segment_length,  # tamanho do segmento de DNA
            need_log=False,  # mostra o log do processo
            genome_path=None  # caminho para o genoma para procurar nullômeros
        )

        with open(dna_output_path, "r") as dna_file:
            data_final.append(dna_file.read())

        # Log para cada parte
        print("Parte {}: Tamanho = {} nucleotídeos".format(idx + 1, len(data_final[-1])))

        # Leitura dos dados para preencher o log e a tabela
        with open(input_path, "r") as original_file, open(dna_output_path, "r") as dna_file:
            original_content = original_file.read()
            dna_content = dna_file.read()
            original_bits = "".join(format(ord(c), "08b") for c in original_content)

        # Escreve no log consolidado
        with open(log_path, "a") as log_file:
            log_file.write("## Parte {}\n\n".format(idx + 1))
            log_file.write("### Informações Gerais\n")
            log_file.write("- Tamanho da sequência (nucleotídeos): {}\n".format(len(dna_content)))
            log_file.write("- Tamanho original (bits): {}\n".format(len(original_bits)))
            log_file.write("- Bits por nucleotídeo: {}\n\n".format(len(original_bits) / len(dna_content)))

            log_file.write("### Detalhamento\n")
            log_file.write("#### Seqüência de Índices\n")
            segment_count = len(dna_content) // segment_length
            for segment_idx in range(segment_count):
                dna_segment = dna_content[segment_idx * segment_length:(segment_idx + 1) * segment_length]
                bit_start = segment_idx * (len(original_bits) // segment_count)
                bit_end = (segment_idx + 1) * (len(original_bits) // segment_count)
                original_segment_bits = original_bits[bit_start:bit_end]

                # Adiciona ao log
                log_file.write("{}. **Índice: {}**\n".format(segment_idx + 1, segment_idx))
                log_file.write("   - DNA: {}\n".format(dna_segment))
                log_file.write("   - Bits originais: {}\n".format(original_segment_bits))

                # Preenche a tabela com as informações do segmento
                index_seq = dna_segment[:16]  # Exemplo: 16 nt para index
                dna_data_seq = dna_segment[16:-16]  # Parte de dados
                full_sequence = "{}{}{}{}".format(primer_5, dna_data_seq, index_seq, primer_3)

                tabela_dados.append({
                    "Indice": segment_idx,
                    "parte": "parte_{}".format(idx + 1),
                    "nt": len(dna_segment),
                    "bits": len(original_segment_bits),
                    "bits/nt": len(original_segment_bits) / len(dna_segment),
                    "primers5": primer_5,
                    "DNA-data": dna_data_seq,
                    "Index": index_seq,
                    "primer3": primer_3,
                    "Sequencia completa": full_sequence,
                    "if original": original_content[bit_start // 8:bit_end // 8],
                    "Bits originais": original_segment_bits
                })

            log_file.write("\n---\n\n")

    # Finaliza a codificação
    print("Log consolidado salvo em {}".format(log_path))

    ### Salvando o arquivo final codificado em DNA
    dna_final_path = os.path.join(dna_dir, "dna_final.dna")
    with open(dna_final_path, "w") as dna_final_file:
        dna_final_file.write("".join(data_final))
    print("DNA final salvo em {}".format(dna_final_path))
    print("Tamanho final em nucleotídeos do Arquivo: {}".format(len("".join(data_final))))

    # Gera o DataFrame para a tabela
    df_tabela = pd.DataFrame(tabela_dados)

    # Salva a tabela em um arquivo Excel
    output_table_path = os.path.join(output_dir, "tabela_segmentos.xlsx")
    df_tabela.to_excel(output_table_path, index=False)

    print("Tabela de segmentos salva em {}".format(output_table_path))


    ## Decodificação / Decoding
    print("DECODIFICAÇÃO DAS PARTES DE DNA EM TEXTO...")
    decoded_final= []
    dna_idx = 0
    for dna_file in sorted(os.listdir(output_partes_dir), key=lambda x: int(x.split("_")[-1].split(".")[0])):
        if dna_file.endswith(".dna"):
            dna_idx += 1
            dna_input_path = os.path.join(output_partes_dir, dna_file)
            model_path = dna_input_path.replace(".dna", ".pkl")
            decoded_output_path = os.path.join(decoded_dir, "decoded_parte_{}.txt".format(dna_idx))

            print("Decodificando Parte {}: {} -> {} usando {}".format(dna_idx, dna_input_path, decoded_output_path, model_path))

            ### Decodificando a sequência de DNA
            pipeline.decode(
            model_path=model_path, #ENTENDI O MODELO É SÓ OS PARAMETROS que foram passados para o YYC, podia ter a informação de index junto
            input_path=dna_input_path,
            output_path=decoded_output_path,
            has_index=True, 
            need_log=False
            # verify=None #méotod de correção de erro deveria vir aqui para fazer sentido
            )

            with open(decoded_output_path, "r") as decoded_file:
                decoded_content = decoded_file.read()
                decoded_final.append(decoded_content)
                print("Parte {} decodificada: {}".format(dna_idx, decoded_content[:50]))  


    ## Junta todas as partes decodificadas em um único arquivo
    ### Primeiro ordena as partes decodificadas para garantir a ordem correta
    decoded_files = sorted(os.listdir(decoded_dir), key=lambda x: int(x.split("_")[-1].split(".")[0]))
    decoded_final_path = os.path.join(output_dir, "decoded_final.txt")
    with open(decoded_final_path, "w") as final_file:
        for decoded_file in decoded_files:
            decoded_file_path = os.path.join(decoded_dir, decoded_file)
            with open(decoded_file_path, "r") as file:
                final_file.write(file.read())

    print("Texto final decodificado salvo em {}".format(decoded_final_path))



    ## Análise de comparação dos arquivos originais e decodificados
    ### Comparação dos arquivos???
    with open(input_file_path, "r") as original_file, open(decoded_final_path, "r") as decoded_file:
        original_text = original_file.read()
        decoded_text = decoded_file.read()
    print("Texto original == Texto decodificado?? ", original_text == decoded_text)

    """
    matrix_1, _ = data_handle.read_binary_from_all(read_file_path, 120, False)
    matrix_2, _ = data_handle.read_binary_from_all(write_file_path, 120, False)
    print("source digital file == target digital file: " + str(matrix_1 == matrix_2))
    ### Printar o tamanho do arquivo original e do arquivo codificado em DNA
    read_file_path = open(read_file_path, "rb").read()
    print("original file size", len(read_file_path))
    dna_path = open(dna_path, "rb").read()
    print("DNA file size", len(dna_path))
    ### printar a taxa de informação
    print("bits per base: ", len(read_file_path)*8 / len(dna_path))
    """