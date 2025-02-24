## **1. Usando a Interface Web do STRINGdb**

A versão web do **STRING** é acessível em https://string-db.org/.

### **1.1 Busca por uma Proteína**

1. Acesse o site e, na barra de pesquisa, insira o nome da proteína, gene ou identificador (ex.: *TP53* para humanos).
2. Selecione a espécie correspondente no menu suspenso.
3. Clique em **Search** para visualizar a rede de interações.

### **1.2 Visualizando Redes de Interação**

•	O resultado principal exibe um **grafo** das interações proteína-proteína com diferentes cores representando fontes de evidência (ex.: interações experimentais, coexpressão, etc.).

•	Passe o cursor sobre os nós para ver informações detalhadas das proteínas.

•	Clique nos nós para explorar conexões mais profundas.

### **1.3 Download de Dados**

•	Você pode baixar as redes de interação no formato **TSV**, **CSV** ou **XGMML** para análises em outros softwares, como Cytoscape.

### **1.4 Análise Funcional**

•	Utilize as guias de **Enriquecimento Funcional** para analisar categorias de GO (*Gene Ontology*), vias de KEGG e associação com doenças.

## **2. Usando o Pacote STRINGdb no R**

Se você tem listas de genes ou proteínas e quer automatizar a análise, pode usar o pacote **STRINGdb** no R.

### **2.1 Instalação do Pacote**

Primeiro, instale e carregue o pacote:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("STRINGdb")
library(STRINGdb)
```

### **2.2 Inicializando a Conexão com STRINGdb**

O primeiro passo é inicializar um objeto STRINGdb com a versão do banco de dados, organismo e a escala de confiança mínima para as interações:

```
string_db <- STRINGdb$new(
    version = "11.5",  # Última versão
    species = 9606,    # Código taxonômico para Homo sapiens (substitua conforme necessário)
    score_threshold = 400,  # Confiança mínima (0 a 1000)
    input_directory = ""
)
```

Para outras espécies, consulte a lista de códigos taxonômicos em:

https://string-db.org/cgi/access?species_text=all

### **2.3 Mapeando Genes para IDs do STRING**

Se você tem uma lista de genes/proteínas, primeiro precisa mapear para os IDs do STRING:

```
# Exemplo de tabela com genes
genes <- data.frame(
    gene = c("TP53", "BRCA1", "BRCA2", "EGFR", "MYC")
)

# Mapeando para os IDs do STRING
mapped_genes <- string_db$map(genes, "gene", removeUnmappedRows = TRUE)
head(mapped_genes)
```

### **2.4 Construindo Redes de Interação**

Com os genes mapeados, podemos gerar uma rede de interações:

```
# Criando um grafo de interações
network <- string_db$get_interactions(mapped_genes$STRING_id)

# Visualizando as primeiras interações
head(network)
```

Se quiser visualizar a rede no próprio STRINGdb:

```
string_db$plot_network(mapped_genes$STRING_id)
```

### **2.5 Análise Funcional (Enriquecimento de GO e KEGG)**

Para obter enriquecimento de termos de **Gene Ontology (GO)** ou **vias de KEGG**:

```
# Análise de enriquecimento GO
enrichment_GO <- string_db$get_enrichment(mapped_genes$STRING_id, category = "Process")
head(enrichment_GO)

# Análise de enriquecimento KEGG
enrichment_KEGG <- string_db$get_enrichment(mapped_genes$STRING_id, category = "KEGG")
head(enrichment_KEGG)
```

### **2.6 Exportando os Resultados**

Para salvar os resultados como arquivos CSV:

```
write.csv(network, "network_interactions.csv", row.names = FALSE)
write.csv(enrichment_GO, "GO_enrichment.csv", row.names = FALSE)
write.csv(enrichment_KEGG, "KEGG_enrichment.csv", row.names = FALSE)
```

## **3. Integração com Cytoscape**

Se quiser visualizar os dados no **Cytoscape**, pode exportar a rede para um formato compatível:

```
write.table(network, "network_for_cytoscape.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
```

No **Cytoscape**, use a opção **File → Import → Network from Table** e carregue o arquivo.

## **Conclusão**

O **STRINGdb** é uma ferramenta poderosa para explorar interações proteína-proteína e analisar funções biológicas associadas. A interface web é útil para visualizações rápidas, enquanto o pacote **STRINGdb no R** permite análises automatizadas e integração com outros dados genômicos.
