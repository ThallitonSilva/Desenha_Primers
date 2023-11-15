import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import primer3

import io


def fasta_genes(fasta):
    stringio = io.StringIO(fasta.getvalue().decode("utf-8"))

    fasta_dos_genes = []
    for record in SeqIO.parse(stringio, 'fasta'):
        fasta_dos_genes.append(SeqRecord(record.seq.upper(),
                                         record.id,
                                         # .split('_cds_')[1].rsplit('_', maxsplit=1)[0],
                                         name='',
                                         description=''))

    return fasta_dos_genes


def cria_listas_genes(lista_genes, fasta):
    stringio = io.StringIO(fasta.getvalue().decode("utf-8"))

    fasta_dos_genes = []
    fasta_total = []

    for record in SeqIO.parse(stringio, "fasta"):

        rec = str(record.id)

        if rec in lista_genes:
            fasta_dos_genes.append(SeqRecord(record.seq.upper(),
                                             record.id,
                                             name='',
                                             description=''))

        fasta_total.append(SeqRecord(record.seq.upper(),
                                     record.id,
                                     name='',
                                     description=''))

    return fasta_dos_genes, fasta_total


def faz_primers(fasta_genes_interesse,
                num_primers,
                prod_min, prod_max,
                p_size_min, p_size_max, p_size_opt,
                p_gc_min, p_gc_max, p_gc_opt,
                p_tm_min, p_tm_max, p_tm_opt):
    result = dict()

    for num, ident in enumerate(fasta_genes_interesse):
        res = primer3.bindings.design_primers(
            {
                'SEQUENCE_ID': fasta_genes_interesse[num].id,
                'SEQUENCE_TEMPLATE': str(fasta_genes_interesse[num].seq),
            },
            {
                'PRIMER_TASK': 'generic',
                'PRIMER_NUM_RETURN': num_primers,
                'PRIMER_PRODUCT_SIZE_RANGE': [prod_min, prod_max],

                'PRIMER_MIN_SIZE': p_size_min,
                'PRIMER_MAX_SIZE': p_size_max,
                'PRIMER_OPT_SIZE': p_size_opt,

                'PRIMER_GC_CLAMP': 1,

                'PRIMER_MIN_GC': p_gc_min,
                'PRIMER_MAX_GC': p_gc_max,
                'PRIMER_OPT_GC_PERCENT': p_gc_opt,

                'PRIMER_MIN_TM': p_tm_min,
                'PRIMER_MAX_TM': p_tm_max,
                'PRIMER_OPT_TM': p_tm_opt,

                'PRIMER_MAX_SELF_ANY': 5,
                'PRIMER_MAX_SELF_END': 3,

                'PRIMER_PAIR_MAX_COMPL_ANY': 5,
                'PRIMER_PAIR_MAX_COMPL_END': 3,
            })
        result[fasta_genes_interesse[num].id] = res

    return result


def organiza_primers(dict_primers):

    lista_primers_info = []

    # Intera sobre cada gene (identificador) e os primers desenhados (result)
    for identificador, result in dict_primers.items():

        # Contagem do numero de primers que foi desenhado para cada gene
        count_primer = result["PRIMER_PAIR_NUM_RETURNED"]

        # Itera sobre cada par de primer
        for num in range(count_primer):
            forward = result[f"PRIMER_LEFT_{num}_SEQUENCE"]
            reverse = result[f"PRIMER_RIGHT_{num}_SEQUENCE"]

            len_f = result[f"PRIMER_LEFT_{num}"][1]
            len_r = result[f"PRIMER_RIGHT_{num}"][1]

            local_anela_f = result[f"PRIMER_LEFT_{num}"][0]
            local_anela_r = result[f"PRIMER_RIGHT_{num}"][0]

            tm_f = result[f"PRIMER_LEFT_{num}_TM"]
            tm_r = result[f"PRIMER_RIGHT_{num}_TM"]

            gc_f = result[f"PRIMER_LEFT_{num}_GC_PERCENT"]
            gc_r = result[f"PRIMER_RIGHT_{num}_GC_PERCENT"]

            product = result[f"PRIMER_PAIR_{num}_PRODUCT_SIZE"]

            lista_primers_info.append([identificador,
                                       forward, reverse,
                                       len_f, len_r,
                                       tm_f, tm_r,
                                       gc_f, gc_r,
                                       local_anela_f, local_anela_r,
                                       product])

    return lista_primers_info


# Crio um DataFrame a partir da lista de primers que passaram no teste
def faz_tabela(dict_primers):
    return pd.DataFrame(dict_primers, columns=['Primer_ID',
                                               'Forward', 'Reverse',
                                               'Tamanho_Forward', 'Tamanho_Reverse',
                                               'TM_Forward', 'TM_Reverse',
                                               'GC_Forward', 'GC_Reverse',
                                               'Anelamento_Forward', 'Anelamento_Reverse',
                                               'Tamanho do fragmento gerado'])


def make_excel(table):
    buffer = io.BytesIO()

    with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:
        table.to_excel(writer, sheet_name=f'Teste_Primers', index=False)

    # writer.save()

    return buffer
