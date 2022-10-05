# hse22_hw1
# 1 часть
1. Создание ссылок
```
ln -s /usr/share/data-minor-bioinf/assembly/oil_R1.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oil_R2.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R1_001.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R2_001.fastq
```
2. Выбор случайных чтений
```
seqtk sample -s724 oil_R1.fastq 5000000 > R1_paired_end.fastq
seqtk sample -s724 oil_R2.fastq 5000000 > R2_paired_end.fastq
seqtk sample -s724 oilMP_S4_L001_R1_001.fastq 1500000 > R1_mate_pairs.fastq
seqtk sample -s724 oilMP_S4_L001_R2_001.fastq 1500000 > R2_mate_pairs.fastq
```
3. Оценка с помощью FastQC
```
mkdir fastqc
ls R*.fastq | xargs -P 4 -tI{} fastqc -o fastqc {}
```
4. Отчет с помощью MultiQC
```
mkdir multiqc
multiqc -o multiqc fastqc

```
5. Обрезание чтений
```
platanus_trim R1_paired_end.fastq R2_paired_end.fastq
platanus_internal_trim R1_mate_pairs.fastq R2_mate_pairs.fastq
```
6. Оценка обрезанных с помощью FastQC
```
ls R* | xargs -P 4 -tI{} fastqc -o fastqc_trimmed {}
```
7. Отчет об обрезанных с помощью MultiQC
```
multiqc -o multiqc_trimmed fastqc_trimmed
```
8. Сбор контиг с помощью “platanus assemble”
```
time platanus assemble -o Poil -f R1_paired_end.fastq.trimmed R2_paired_end.fastq.trimmed 2> assemble.log
```
9. Сбор скаффолдов
```
time platanus scaffold -o Poil -c Poil_contig.fa -IP1 R1_paired_end.fastq.trimmed R2_paired_end.fastq.trimmed -OP2 R1_mate_pairs.fastq.int_trimmed R2_mate_pairs.fastq.int_trimmed 2> scaffold.log
```
10. Уменьшение числа промежутков
```
time platanus gap_close -o Poil -c Poil_scaffold.fa -IP1 R1_paired_end.fastq.trimmed R2_paired_end.fastq.trimmed -OP2 R1_mate_pairs.fastq.int_trimmed R2_mate_pairs.fastq.int_trimmed 2> gapclose.log
```
## Отчеты multiQC
Для исходных чтений
![avatar](/images/general_1.png)
![avatar](/images/per_sequence_1.png)
Для подрезанных
![avatar](/images/general_2.png)
![avatar](/images/per_sequence_2.png)
## 2 часть
1. Импорт библиотек
```
import re
```
2. Функция для получения данных
```
def analysis(path):
    file = open(path, "r")
    lengths = []
    total_length, quantity, max_length, length, score = 0, 0, 0, 0, 0
    maximum_sequence, current_sequence = "", ""
    for line in file:
        if (line[0] == ">"):
            if quantity != 0:
                lengths.append(length)
            quantity += 1
            if length >= max_length:
                max_length = length
                maximum_sequence = current_sequence
            current_sequence = ""
            length = 0
        else:
            current_sequence += line.strip()
            length += len(line.strip())
            total_length += len(line.strip())
     
    lengths.sort(reverse = True) 
    for i in lengths:
        score += i
        if score >= total_length / 2:
          print(f"Общее количество: {quantity},\nОбщая длина: {total_length},\n\
Максимальная длина: {max_length},\nN50: {i}\n")
          break
    return maximum_sequence
```
3. Контиги
```
contig = analysis("Poil_contig.fa")
```
Общее количество: 606,
Общая длина: 3923405,
Максимальная длина: 179304,
N50: 48054

4. Скаффолды
```
scaffolds = analysis("Poil_scaffold.fa")
```
Общее количество: 71,
Общая длина: 3871958,
Максимальная длина: 3831358,
N50: 3831358

5. Подсчет гэпов для необрезанных чтений
```
print("Гэпы:")
print(f"Общая длина: {scaffolds.count('N')}")
scaffolds = re.sub(r"N{2,}", "N", scaffolds)
print(f"Количество: {scaffolds.count('N')}")
```
Гэпы:
Общая длина: 6338
Количество: 64

6. Подсчет гэпов для обрезанных чтений
```
scaffold_closed = analysis("Poil_gapClosed.fa")
print(f"Общая длина для обрезанных чтений: {scaffold_closed.count('N')}")
scaffold_closed = re.sub(r"N{2,}", "N", scaffold_closed)
print(f"Количество для обрезанных чтений: {scaffold_closed.count('N')}")
```
Общая длина для обрезанных чтений: 1740
Количество для обрезанных чтений: 8