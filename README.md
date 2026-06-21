# zuzalysin-homologues-database

# Baza homologów zuzalizyny — opis kolumn

Baza zawiera homologi i paralogi **zuzalizyny** (białko Q7MTD8 z *Porphyromonas gingivalis*)
znalezione w genomach bakterii. **Jeden wiersz = jeden paralog** w danym genomie.
Paralogi z tego samego genomu mają identyczne wartości w kolumnach „genomowych”
(taksonomia, statystyki assembly, CheckM, T9SS) i różne w kolumnach „paralogowych”
(sekwencja, domeny, peptyd sygnałowy).

Genomy grupowane są po kolumnie **`Assembly`** (numer GCF).

---

## Trafienie BLAST i taksonomia

| Kolumna | Co przechowuje |
|---|---|
| `Scientific Name` | Nazwa gatunku (organizmu), z którego pochodzi trafienie. |
| `Taxid` | Identyfikator taksonomiczny NCBI (Taxonomy ID). |
| `Max Score` | Najwyższy bit-score pojedynczego dopasowania (HSP) w BLAST. |
| `Total Score` | Sumaryczny bit-score wszystkich dopasowań danego trafienia. |
| `Query Cover` | % długości sekwencji zapytania (zuzalizyny) pokryty dopasowaniem. |
| `E value` | Wartość E — statystyczna istotność dopasowania BLAST (im mniejsza, tym lepiej). |
| `Per. ident` | % identyczności aminokwasów w dopasowaniu do zuzalizyny. |
| `Acc. Len` | Długość rekordu (sekwencji nukleotydowej) z trafienia. |
| `Accession_ID` | Numer dostępu rekordu nukleotydowego w NCBI (np. `NC_010729`). |
| `Accession_URL` | Bezpośredni link do rekordu w NCBI. |
| `paralogi_blast` | Liczba paralogów znalezionych przez BLAST w tym genomie. |
| `szacowana_liczba_paralogów` | Szacowana liczba paralogów w genomie. |

## Genom / assembly (NCBI Datasets)

| Kolumna | Co przechowuje |
|---|---|
| `Assembly` | Numer złożenia genomu (GCF…). **Klucz grupujący paralogi w genom.** |
| `Assembly Name` | Nazwa złożenia (np. `ASM1050v1`). |
| `ANI Check status` | Wynik weryfikacji tożsamości genomu metodą ANI (np. `OK`). |
| `Organism Infraspecific Names Strain` | Nazwa/oznaczenie szczepu. |
| `Annotation Pipeline` | Pipeline anotacji genomu (np. NCBI PGAP). |
| `Annotation Release Date` | Data anotacji genomu. |
| `Assembly Stats Total Sequence Length` | Całkowita długość genomu (w parach zasad, bp). |
| `Assembly Stats Total Number of Chromosomes` | Liczba chromosomów. |
| `Assembly Stats Number of Contigs` | Liczba kontigów. |
| `Assembly Level` | Poziom kompletności złożenia (Complete Genome / Chromosome / Scaffold / Contig). |
| `Assembly Stats Contig N50` | N50 kontigów (miara ciągłości złożenia). |
| `Assembly Stats Scaffold N50` | N50 scaffoldów. |
| `Assembly Stats Number of Scaffolds` | Liczba scaffoldów. |
| `Assembly Stats GC Percent` | Zawartość par GC w genomie (%). |
| `Assembly Sequencing Tech` | Technologia sekwencjonowania (np. PacBio). |
| `Assembly BioProject Accession` | Numer projektu BioProject. |
| `Assembly BioSample Accession` | Numer próbki BioSample. |
| `Annotation Count Gene Total` | Łączna liczba genów w genomie. |
| `Annotation Count Gene Protein-coding` | Liczba genów kodujących białka. |
| `Annotation Count Gene Pseudogene` | Liczba pseudogenów. |
| `Type Material Display Text` | Informacja, czy genom pochodzi z materiału typowego szczepu. |
| `CheckM marker set` | Zestaw markerów użyty przez CheckM do oceny jakości. |
| `CheckM completeness` | Szacowana kompletność genomu wg CheckM (%). |
| `CheckM contamination` | Szacowana kontaminacja genomu wg CheckM (%). |

## System sekrecji typu IX (T9SS)

| Kolumna | Co przechowuje |
|---|---|
| `SprA_present` / `PorN_present` / `PorU_present` / `PorV_present` | Czy w genomie wykryto dany komponent T9SS (`TAK`/`NIE`). |
| `SprA_evalue` / `PorN_evalue` / `PorU_evalue` / `PorV_evalue` | e-value dopasowania danego komponentu. |
| `SprA_count` / `PorN_count` / `PorU_count` / `PorV_count` | Liczba znalezionych kopii danego komponentu. |
| `T9SS_components_found` | Liczba znalezionych komponentów T9SS (0–4). |
| `T9SS_status` | Status systemu T9SS: `complete`, `partial`, `partial_no_translocon`, `absent`. |

## Paralog — domeny i motyw katalityczny

| Kolumna | Co przechowuje |
|---|---|
| `Paralog_ID` | Identyfikator paralogu w formacie `paralog_<nr>_pos<start>_<end>` (pozycje w genomie). |
| `czy_jest_HEXXHXXGXXH` | Czy obecny jest motyw katalityczny metaloproteazy `HExxHxxGxxH` (`TAK`/`NIE`). |
| `DUF4953_evalue` / `DUF4953_score` | e-value i score dopasowania domeny **DUF4953** (model HMM). |
| `DUF5117_evalue` / `DUF5117_score` | e-value i score dopasowania domeny **DUF5117** (domena peptydazy). |
| `DUF5118_evalue` / `DUF5118_score` | e-value i score dopasowania domeny **DUF5118**. |
| `czy_jest_DUF4953` / `czy_jest_DUF5117` / `czy_jest_DUF5118` | Czy dana domena została wykryta (`TAK`/`NIE`). |

## Paralog — peptyd sygnałowy (SignalP)

| Kolumna | Co przechowuje |
|---|---|
| `Prediction` | Typ peptydu sygnałowego przewidziany przez SignalP: `OTHER` (brak), `SP`, `LIPO`, `TAT`, `TATLIPO`, `PILIN`. |
| `OTHER` | Prawdopodobieństwo braku peptydu sygnałowego. |
| `SP(Sec/SPI)` | Prawdopodobieństwo standardowego peptydu sygnałowego (szlak Sec, peptydaza SPI). |
| `LIPO(Sec/SPII)` | Prawdopodobieństwo peptydu lipoproteinowego (Sec, SPII). |
| `TAT(Tat/SPI)` | Prawdopodobieństwo peptydu szlaku Tat (SPI). |
| `TATLIPO(Tat/SPII)` | Prawdopodobieństwo peptydu Tat lipoproteinowego (SPII). |
| `PILIN(Sec/SPIII)` | Prawdopodobieństwo peptydu pilin (SPIII). |
| `CS Position` | Pozycja miejsca cięcia peptydu sygnałowego (np. `CS pos: 23-24. Pr: 0.95`). Odnosi się do **pełnej** sekwencji. |

## Paralog — sekwencje i identyfikator białka

| Kolumna | Co przechowuje |
|---|---|
| `sequence` | Sekwencja aminokwasowa **fragmentu dopasowania BLAST** — może nie zawierać N-końca, w tym peptydu sygnałowego. |
| `sekwencje_cale` | **Pełna sekwencja aminokwasowa białka** (z N-końcem i peptydem sygnałowym). Podstawa wyświetlania w aplikacji oraz analizy SignalP. |
| `ID_bialka` | Numer dostępu białka w NCBI (np. `WP_012458658.1`). |

---

## Jak aplikacja wyświetla sekwencję

1. Jeśli `sekwencje_cale` ma wartość → pokazywana jest **pełna sekwencja**, z zaznaczonym
   **peptydem sygnałowym** (wg `CS Position`) oraz **motywem `HExxHxxGxxH`**.
2. Jeśli `sekwencje_cale` jest pusta → pokazywana jest `sequence` (sam fragment dopasowania),
   pojawia się komunikat **„UWAGA: to nie cała sekwencja, tylko dopasowanie”**, a peptyd
   sygnałowy **nie jest zaznaczany** (jego pozycja dotyczy pełnego białka, nie fragmentu).

## Cechy wyliczane przez aplikację (nie ma ich w pliku CSV)

- `architecture` — klasa budowy paralogu (od „Pełna (z peptydem sygnałowym)” do „Bez motywu katalitycznego”).
- `completeness_score` (0–5) — suma: motyw + DUF4953 + DUF5117 + DUF5118 + peptyd sygnałowy.
- `length_aa` — długość sekwencji wyświetlanej (pełnej, gdy jest; inaczej dopasowania).
- `n_DUF` — liczba wykrytych domen DUF (0–3).

# Instrukcja uruchomienia GUI bazy danych

Ta instrukcja prowadzi krok po kroku przez proces uruchomienia graficznego interfejsu użytkownika (GUI) na komputerze, na którym nie jest jeszcze zainstalowane żadne oprogramowanie programistyczne.

---

## Krok 1: Pobierz repozytorium z GitHub

1. Wejdź na stronę repozytorium na GitHubie.
2. Kliknij zielony przycisk **Code**, a następnie wybierz **Download ZIP**.
3. Rozpakuj pobrany plik `.zip` w wygodnym miejscu na komputerze, np. na Pulpicie.

Po rozpakowaniu folder powinien zawierać m.in.:

- `gui.py` — kod aplikacji
- `database.csv` — baza danych
- `tree.svg` — drzewo filogenetyczne
- `pyproject.toml` — lista bibliotek wymaganych przez aplikację
- `uv.lock` — zablokowane wersje bibliotek
- `README.md`

---

## Krok 2: Zainstaluj `uv`

`uv` to narzędzie, które automatycznie zainstaluje odpowiednią wersję Pythona oraz wszystkie potrzebne biblioteki — nie trzeba niczego instalować ręcznie wcześniej.

### Windows

Otwórz **PowerShell** (wpisz w wyszukiwarce systemowej "PowerShell") i wklej:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

### macOS

Otwórz aplikację **Terminal** (Launchpad → Inne → Terminal) i wklej:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Po instalacji **zamknij i ponownie otwórz** okno terminala/PowerShell, aby zmiany zaczęły działać.

Sprawdź, czy instalacja się powiodła:

```bash
uv --version
```

Powinna pojawić się informacja o numerze wersji (np. `uv 0.5.x`).

---

## Krok 3: Przejdź do folderu z aplikacją

W oknie terminala/PowerShell przejdź do folderu, w którym rozpakowano repozytorium. Przykład (jeśli folder nazywa się `baza_zuzalizyny` i znajduje się na Pulpicie):

**Windows:**
```powershell
cd Desktop\baza_zuzalizyny
```

**macOS:**
```bash
cd Desktop/baza_zuzalizyny
```

---

## Krok 4: Zainstaluj zależności

W folderze repozytorium znajdują się dwa pliki istotne dla `uv`:

- **`pyproject.toml`** — lista bibliotek wymaganych przez aplikację oraz wersja Pythona
- **`uv.lock`** — dokładne, zablokowane wersje wszystkich bibliotek, gwarantujące identyczne działanie aplikacji niezależnie od komputera

Dzięki nim wystarczy jedna komenda, aby `uv` automatycznie:
- pobrał odpowiednią wersję Pythona (jeśli nie jest jeszcze zainstalowana),
- utworzył izolowane środowisko,
- zainstalował wszystkie wymagane biblioteki w dokładnie tych wersjach, które zostały przetestowane.

Wpisz w terminalu:

```bash
uv sync
```

Poczekaj, aż proces się zakończy (może to potrwać od kilkudziesięciu sekund do kilku minut, w zależności od szybkości łącza internetowego).

---

## Krok 5: Uruchom aplikację

Po zakończeniu instalacji uruchom aplikację komendą:

```bash
uv run streamlit run gui.py
```

Aplikacja powinna **automatycznie otworzyć się w przeglądarce internetowej** pod adresem podobnym do:

```
http://localhost:8501
```

Jeśli przeglądarka nie otworzy się sama, skopiuj ten adres z terminala i wklej go ręcznie w pasek przeglądarki.

---

## Krok 6: Korzystanie z aplikacji

W uruchomionej aplikacji można:

- przeglądać i filtrować bazę danych organizmów według wybranych kryteriów (np. obecność peptydu sygnałowego, status T9SS),
- przeglądać szczegółowe informacje o poszczególnych paralogach (sekwencja, domeny, motyw katalityczny),
- przeglądać interaktywne drzewo filogenetyczne,
- pobierać wybrane sekwencje w formacie `.fasta`.

---

## Zamykanie aplikacji

Aby zatrzymać działanie aplikacji, wróć do okna terminala/PowerShell i naciśnij jednocześnie klawisze:

```
Ctrl + C
```

Następnie można zamknąć okno terminala.

---

## Najczęstsze problemy

**„uv nie jest rozpoznawany jako polecenie”**
Terminal/PowerShell został otwarty przed zainstalowaniem `uv` lub przed jego dodaniem do ścieżki systemowej. Zamknij i otwórz terminal ponownie, a następnie spróbuj jeszcze raz.

**Polecenie `cd` nie znajduje folderu**
Upewnij się, że ścieżka do folderu jest poprawna. W razie wątpliwości można przeciągnąć folder z Eksploratora plików / Findera bezpośrednio do okna terminala — ścieżka wstawi się automatycznie.

**Strona w przeglądarce się nie otwiera**
Sprawdź w oknie terminala, czy pojawił się adres zaczynający się od `http://localhost`, i wklej go ręcznie do przeglądarki.

**Po zamknięciu terminala aplikacja przestaje działać**
To normalne — aplikacja działa tylko, gdy terminal jest otwarty. Aby ją ponownie uruchomić, wystarczy powtórzyć Krok 5 (nie trzeba powtarzać Kroku 4).


# Licencja

Ten projekt jest udostępniony na licencji 
[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).

Oznacza to, że można:
- swobodnie korzystać z danych i kodu w celach niekomercyjnych
- modyfikować i rozwijać projekt

pod warunkiem podania autora (cytowania) i niewykorzystywania 
do celów komercyjnych.
