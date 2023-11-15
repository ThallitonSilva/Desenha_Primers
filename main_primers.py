import streamlit as st
import streamlit_ext as ste
from funcs_primers import *

#################### ----------------------- ###################################### ----------------------- ############

st.set_page_config(page_title='Desenho dos Primers', layout='wide')

hide_menu_style = """
        <style>
        #MainMenu {visibility: hidden; }
        footer {visibility: hidden;}
        </style>
        """

st.markdown(hide_menu_style, unsafe_allow_html=True)

#################### ----------------------- ###################################### ----------------------- ############

title_col1, title_col2, title_col3, title_col4 = st.columns((1, 2, 2, 1))
col1, col2, col3, col4 = st.columns((1, 2, 2, 1))

title_col2.title('Desenho dos primers')

primers = col2.file_uploader('Insira um arquivo excel com os ids dos genes que se deseja desenhar os primers',
                             type='xlsx',
                             accept_multiple_files=False)

genome = col3.file_uploader('Insira o genoma (formato FASTA) em que os primers serão baseados',
                            type='fasta',
                            accept_multiple_files=False)

exp = col2.expander("""A tabela deve, obrigatoriamente, 
ter a coluna Gene_ID, com o mesmo ID do gene presente no genoma inserido. Exemplo:""")

exemplo = [['Gene_ID'],
           ['AT1G60590'],
           ['AT4G37990'],
           ['AT2G33380']]

ex = pd.DataFrame(exemplo[1:], columns=exemplo[0])
exp.write(ex)

if primers and genome:

    table_primers = pd.read_excel(primers)

    col2.write(table_primers)

    col2.write(f'Você deseja desenhar primers para {table_primers.shape[0]} genes')

    list_fasta = fasta_genes(genome)

    col3.write(f'O arquivo fasta inserido possui {len(list_fasta)} sequencias')

    if 'Gene_ID' not in table_primers.columns:
        st.write('O nome das colunas estão errados. Por favor, arrume!')
        st.write('Veja o exemplo fornecido.')
        col2.write(table_primers)

    else:
        sidebar = st.sidebar

        sidebar.write('Parameters')

        sidebar.divider()

        sidebar.write('Number')

        num_primers = sidebar.slider('Select primer pair numbers to draw',
                                     1, 200,
                                     10)

        sidebar.divider()

        sidebar.write('Length of products')

        prod_min, prod_max = sidebar.slider('Select the min/max length of the products of PCR',
                                            75, 500,
                                            (80, 150))

        sidebar.divider()

        sidebar.write('Length of primers')

        p_size_min, p_size_max = sidebar.slider('Select the min/max length of primers',
                                                15, 35,
                                                (18, 22))

        p_size_opt = sidebar.select_slider('Select the optimal lenght',
                                           options=range(p_size_min, p_size_max+1),
                                           value=int(p_size_min + ((p_size_max - p_size_min) / 2)))

        sidebar.divider()

        sidebar.write('GC content')

        p_gc_min, p_gc_max = sidebar.slider('Select the min/max %GC for the primers',
                                            0, 100,
                                            (40, 60))

        p_gc_opt = sidebar.select_slider('Select the optimal %GC',
                                         options=range(p_gc_min, p_gc_max+1),
                                         value=int(p_gc_min + ((p_gc_max - p_gc_min) / 2)))

        sidebar.divider()

        sidebar.write('Temperature')

        p_tm_min, p_tm_max = sidebar.slider('Select the min/max TM of the primers',
                                            30, 90,
                                            (50, 70))

        p_tm_opt = sidebar.select_slider('Select the optimal TM',
                                         options=range(p_tm_min, p_tm_max+1),
                                         value=int(p_tm_min + ((p_tm_max - p_tm_min) / 2)))

        sidebar.divider()

        run_app = sidebar.button('Run')


    if run_app:

        st.write('Iniciando o Desenho dos primers...')

        ids = table_primers['Gene_ID'].to_list()

        fasta_genes, fasta_tot = cria_listas_genes(ids, genome)

        st.write('Fastas gerados...')

        primers = faz_primers(fasta_genes,
                              num_primers,
                              prod_min, prod_max,
                              p_size_min, p_size_max, p_size_opt,
                              p_gc_min, p_gc_max, p_gc_opt,
                              p_tm_min, p_tm_max, p_tm_opt)

        st.write('Primers desenhados...')

        print(primers)

        organizado = organiza_primers(primers)

        tabela = faz_tabela(organizado)

        st.write('Tabela Feita...')
        st.write(tabela)

        excel = make_excel(tabela)

        ste.download_button(label="Download Resultados",
                            data=excel,
                            file_name="Desenho_dos_Primers.xlsx",
                            mime="application/vnd.ms-excel")
            
        st.write('Agora vá até: https://teste-primers.streamlit.app/ e coloque essa tabela gerada junto com o genoma para testar os primers desenhados!!')    

#################### ----------------------- ###################################### ----------------------- ############
