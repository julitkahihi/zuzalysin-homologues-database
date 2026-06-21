from __future__ import annotations
import io
import re
import pandas as pd
import numpy as np
import streamlit as st
import streamlit.components.v1 as components
import plotly.express as px
import plotly.graph_objects as go


st.set_page_config(page_title="Baza homologów zuzalizyny",
                   layout="wide", initial_sidebar_state="expanded")

PER_PARALOG_COLS = [
    "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)",
    "TATLIPO(Tat/SPII)", "PILIN(Sec/SPIII)", "CS Position", "Paralog_ID", "sequence",
    "sekwencje_cale", "ID_bialka",
    "czy_jest_HEXXHXXGXXH", "DUF4953_evalue", "DUF4953_score", "DUF5118_evalue",
    "DUF5118_score", "DUF5117_evalue", "DUF5117_score", "czy_jest_DUF4953",
    "czy_jest_DUF5117", "czy_jest_DUF5118",
]

NUMERIC_COLS = [
    "Taxid", "Max Score", "Total Score", "Per. ident", "Acc. Len", "paralogi_blast",
    "szacowana_liczba_paralogów", "Assembly Stats Total Sequence Length",
    "Assembly Stats Total Number of Chromosomes", "Assembly Stats Number of Contigs",
    "Assembly Stats Contig N50", "Assembly Stats Scaffold N50",
    "Assembly Stats Number of Scaffolds", "Assembly Stats GC Percent",
    "Annotation Count Gene Total", "Annotation Count Gene Protein-coding",
    "Annotation Count Gene Pseudogene", "CheckM completeness", "CheckM contamination",
    "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)", "TATLIPO(Tat/SPII)",
    "PILIN(Sec/SPIII)", "SprA_evalue", "SprA_count", "PorN_evalue", "PorN_count",
    "PorU_evalue", "PorU_count", "PorV_evalue", "PorV_count", "T9SS_components_found",
    "DUF4953_evalue", "DUF4953_score", "DUF5118_evalue", "DUF5118_score",
    "DUF5117_evalue", "DUF5117_score", "E value",
]

MOTIF_REGEX = re.compile(r"HE..H..G..H")  # HExxHxxGxxH 

T9SS_ORDER = ["complete", "partial", "partial_no_translocon", "absent"]
T9SS_COLORS = {"complete": "#1f7a6f", "partial": "#d9a441",
               "partial_no_translocon": "#d97b3f", "absent": "#9aa3a0"}
ARCH_ORDER = ["Pełna (z peptydem sygnałowym)", "Pełna domenowo (bez peptydu)",
              "Z motywem, niepełne domeny", "Bez motywu katalitycznego"]
ARCH_COLORS = {ARCH_ORDER[0]: "#1f7a6f", ARCH_ORDER[1]: "#5b8def",
               ARCH_ORDER[2]: "#d9a441", ARCH_ORDER[3]: "#c2503a"}
PALETTE = ["#1f7a6f", "#d98c3f", "#5b8def", "#9b5bbb", "#c2503a", "#3fb6a8",
           "#e0b94d", "#6b7280", "#8a9a5b", "#b5651d"]
BOOL_COLS = ["czy_jest_HEXXHXXGXXH", "czy_jest_DUF4953", "czy_jest_DUF5117",
             "czy_jest_DUF5118"]

PLOTLY_TEMPLATE = "plotly_white"



def _read_csv(data) -> pd.DataFrame:

    df = pd.read_csv(data, dtype=str, skipinitialspace=True, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
    df = df.replace({"": np.nan, "nan": np.nan, "NaN": np.nan, "None": np.nan})
    return df


def _coerce(df: pd.DataFrame) -> pd.DataFrame:
    for c in NUMERIC_COLS:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if "Query Cover" in df.columns:
        df["Query Cover (%)"] = pd.to_numeric(
            df["Query Cover"].astype(str).str.replace("%", "", regex=False),
            errors="coerce")
    return df


def _bool(df: pd.DataFrame, col: str) -> pd.Series:
    if col in df.columns:
        return df[col].map({"TAK": True, "NIE": False})
    return pd.Series([np.nan] * len(df), index=df.index)


def _derive_paralog_features(df: pd.DataFrame) -> pd.DataFrame:
    df["has_HEXXH"] = _bool(df, "czy_jest_HEXXHXXGXXH").fillna(False)
    df["has_DUF4953"] = _bool(df, "czy_jest_DUF4953").fillna(False)
    df["has_DUF5117"] = _bool(df, "czy_jest_DUF5117").fillna(False)
    df["has_DUF5118"] = _bool(df, "czy_jest_DUF5118").fillna(False)
    df["n_DUF"] = (df["has_DUF4953"].astype(int) + df["has_DUF5117"].astype(int)
                   + df["has_DUF5118"].astype(int))
    df["has_all_DUF"] = df["has_DUF4953"] & df["has_DUF5117"] & df["has_DUF5118"]
    if "Prediction" in df.columns:
        df["has_signal"] = df["Prediction"].notna() & (df["Prediction"] != "OTHER")
    else:
        df["has_signal"] = False
    df["completeness_score"] = (df["has_HEXXH"].astype(int)
                                + df["has_DUF4953"].astype(int)
                                + df["has_DUF5117"].astype(int)
                                + df["has_DUF5118"].astype(int)
                                + df["has_signal"].astype(int))

    def classify(r):
        if r.has_HEXXH and r.has_all_DUF and r.has_signal:
            return ARCH_ORDER[0]
        if r.has_HEXXH and r.has_all_DUF:
            return ARCH_ORDER[1]
        if r.has_HEXXH:
            return ARCH_ORDER[2]
        return ARCH_ORDER[3]
    df["architecture"] = df.apply(classify, axis=1)

    # Sekwencja do wyświetlania: pełna ("sekwencje_cale") gdy dostępna,
    # inaczej samo dopasowanie BLAST ("sequence").
    if "sekwencje_cale" in df.columns:
        df["has_full_seq"] = (df["sekwencje_cale"].notna()
                              & (df["sekwencje_cale"].astype(str).str.strip() != ""))
    else:
        df["has_full_seq"] = pd.Series(False, index=df.index)
    seq_part = (df["sequence"] if "sequence" in df.columns
                else pd.Series([np.nan] * len(df), index=df.index))
    seq_full = (df["sekwencje_cale"] if "sekwencje_cale" in df.columns
                else pd.Series([np.nan] * len(df), index=df.index))
    df["display_seq"] = seq_full.where(df["has_full_seq"], seq_part)
    # Długość liczona z sekwencji wyświetlanej (pełnej, gdy jest; inaczej dopasowania).
    df["length_aa"] = df["display_seq"].fillna("").astype(str).str.len()

    if "Paralog_ID" in df.columns:
        ext = df["Paralog_ID"].astype(str).str.extract(
            r"paralog_(?P<num>\d+)_pos(?P<start>\d+)_(?P<end>\d+)")
        df["paralog_num"] = pd.to_numeric(ext["num"], errors="coerce")
        df["genome_start"] = pd.to_numeric(ext["start"], errors="coerce")
        df["genome_end"] = pd.to_numeric(ext["end"], errors="coerce")
    return df


def _make_genome_key(df: pd.DataFrame) -> pd.DataFrame:
    key = df["Assembly"] if "Assembly" in df.columns else pd.Series([np.nan] * len(df))
    key = key.copy()
    if "Accession_ID" in df.columns:
        key = key.fillna(df["Accession_ID"])
    if "Scientific Name" in df.columns:
        key = key.fillna(df["Scientific Name"])
    key = key.fillna(pd.Series([f"genome_{i}" for i in range(len(df))], index=df.index))
    df["__genome"] = key
    if "Scientific Name" in df.columns:
        df["Genus"] = df["Scientific Name"].astype(str).str.split().str[0]
    else:
        df["Genus"] = "—"
    sci = df["Scientific Name"] if "Scientific Name" in df.columns else df["__genome"]
    strain = (df["Organism Infraspecific Names Strain"]
              if "Organism Infraspecific Names Strain" in df.columns
              else pd.Series([np.nan] * len(df), index=df.index))
    df["genome_label"] = np.where(strain.notna(), sci + " — " + strain.astype(str), sci)
    return df


def build_genome_table(df: pd.DataFrame) -> pd.DataFrame:
    """Zwija paralogi do poziomu genomu (jeden wiersz na genom)."""
    genome_cols = [c for c in df.columns
                   if c not in PER_PARALOG_COLS
                   and c not in ("has_HEXXH", "has_DUF4953", "has_DUF5117",
                                 "has_DUF5118", "n_DUF", "has_all_DUF", "has_signal",
                                 "completeness_score", "architecture", "length_aa",
                                 "paralog_num", "genome_start", "genome_end",
                                 "has_full_seq", "display_seq")]
    base = df.drop_duplicates("__genome")[genome_cols].set_index("__genome")
    agg = df.groupby("__genome").agg(
        n_paralogs=("Paralog_ID", "size"),
        n_HEXXH=("has_HEXXH", "sum"),
        n_full=("has_all_DUF", "sum"),
        n_secreted=("has_signal", "sum"),
        n_canonical=("architecture", lambda s: int((s == ARCH_ORDER[0]).sum())),
        max_DUF5117=("DUF5117_score", "max"),
        max_DUF4953=("DUF4953_score", "max"),
        mean_len=("length_aa", "mean"),
    )
    g = base.join(agg).reset_index()
    g["pct_HEXXH"] = np.where(g["n_paralogs"] > 0, 100 * g["n_HEXXH"] / g["n_paralogs"], 0)
    return g


@st.cache_data(show_spinner="Przetwarzam dane…")
def process(raw_bytes: bytes):
    df = _read_csv(io.BytesIO(raw_bytes))
    df = _coerce(df)
    df = _make_genome_key(df)
    df = _derive_paralog_features(df)
    genomes = build_genome_table(df)
    return df, genomes


def yn(val) -> str:
    return {True: "✅", False: "❌"}.get(bool(val), "—") if pd.notna(val) else "—"


def fmt(v, nd=1):
    if pd.isna(v):
        return "—"
    if isinstance(v, (int, np.integer)):
        return f"{int(v):,}".replace(",", " ")
    return f"{v:,.{nd}f}".replace(",", " ")


def highlight_sequence_html(seq: str, cs_first: int | None) -> str:
    seq = seq or ""
    n = len(seq)
    flags = [""] * n
    if cs_first:
        for i in range(min(int(cs_first), n)):
            flags[i] = "sig"
    for m in MOTIF_REGEX.finditer(seq):
        for i in range(m.start(), m.end()):
            flags[i] = "mot"
    width = 60
    lines = []
    for i in range(0, n, width):
        chunk = []
        for j in range(i, min(i + width, n)):
            cls = flags[j]
            ch = seq[j]
            chunk.append(f"<span class='{cls}'>{ch}</span>" if cls else ch)
        lines.append(
            f"<div class='seqline'><span class='ln'>{i + 1}</span>"
            f"<span class='res'>{''.join(chunk)}</span></div>")
    return f"<div class='seqbox'>{''.join(lines)}</div>"


def parse_cs_first(cs: str):
    if not isinstance(cs, str):
        return None
    m = re.search(r"CS pos:\s*(\d+)-(\d+)", cs)
    return int(m.group(1)) if m else None


def to_fasta(rows: pd.DataFrame) -> str:
    lines = []
    for _, r in rows.iterrows():
        sci = r.get("Scientific Name", "")
        prot = r.get("ID_bialka", "")
        prot_tag = f"|{prot}" if isinstance(prot, str) and prot.strip() else ""
        header = f">{r['__genome']}|{r.get('Paralog_ID','')}{prot_tag} {sci}".strip()
        # eksportujemy pełną sekwencję, gdy jest dostępna; inaczej dopasowanie
        seq = r.get("display_seq")
        if not isinstance(seq, str) or not seq:
            seq = r.get("sequence") or ""
        lines.append(header)
        for i in range(0, len(seq), 60):
            lines.append(seq[i:i + 60])
    return "\n".join(lines)


def styled_bool_df(df: pd.DataFrame, bool_like_cols):
    out = df.copy()
    for c in bool_like_cols:
        if c in out.columns:
            out[c] = out[c].map({"TAK": "✅", "NIE": "❌", True: "✅", False: "❌"})
    return out


CSS = """
<style>
:root { --pri:#1f7a6f; --ink:#1c2b2a; }
.block-container { padding-top: 1.6rem; max-width: 1500px; }
h1,h2,h3 { color: var(--ink); letter-spacing:-.01em; }
.kpi {
  background: linear-gradient(160deg,#ffffff 0%,#f3f7f6 100%);
  border:1px solid #e3eae8; border-left:4px solid var(--pri);
  border-radius:14px; padding:14px 16px;
  min-height:104px; margin-bottom:14px; box-sizing:border-box;
}
.kpi .lab { font-size:.74rem; text-transform:uppercase; letter-spacing:.06em;
  color:#5d6b69; font-weight:600; }
.kpi .val { font-size:1.7rem; font-weight:700; color:var(--ink); line-height:1.15; }
.kpi .sub { font-size:.78rem; color:#7a8785; }
.badge { display:inline-block; padding:3px 12px; border-radius:999px;
  font-size:.82rem; font-weight:600; color:#fff; }
.detailcard { background:#fbfdfc; border:1px solid #e3eae8; border-radius:14px;
  padding:18px 20px; }
.seqbox { background:#0f1c1a; border-radius:10px; padding:14px 16px;
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
  font-size:.82rem; line-height:1.65; overflow-x:auto; }
.seqbox .seqline { white-space:nowrap; }
.seqbox .ln { display:inline-block; width:3.5em; text-align:right;
  margin-right:1.1em; color:#5b736e; user-select:none; }
.seqbox .res { color:#cfe3df; }
.seqbox .mot { background:#c2503a; color:#fff; border-radius:3px; font-weight:700; }
.seqbox .sig { background:#345; color:#bcd; border-bottom:2px solid #5b8def; }
.small { font-size:.82rem; color:#5d6b69; }
hr { margin:.8rem 0; border:none; border-top:1px solid #e3eae8; }
</style>
"""


def kpi(col, label, value, sub=""):
    col.markdown(
        f"<div class='kpi'><div class='lab'>{label}</div>"
        f"<div class='val'>{value}</div><div class='sub'>{sub}</div></div>",
        unsafe_allow_html=True)


# Podgląd drzewa filogenetycznego z SVG
def pokaz_drzewo_svg(svg_text, height=700):
    """Osadza gotowy SVG z własnym, lekkim zoomem i przesuwaniem napisanym
    w czystym JavaScript — bez żadnych zewnętrznych bibliotek ani linków."""
    i = svg_text.find("<svg")
    if i > 0:
        svg_text = svg_text[i:]
    btn = ("border:1px solid #cdd6d4;background:#fff;cursor:pointer;border-radius:8px;"
           "height:34px;min-width:34px;font-size:16px;")
    html = f"""
    <div style="position:relative;border:1px solid #e3eae8;border-radius:12px;
                overflow:hidden;background:#fff;height:{height}px;">
      <div style="position:absolute;top:10px;right:10px;z-index:5;display:flex;gap:6px;">
        <button id="z_in"  style="{btn}">+</button>
        <button id="z_out" style="{btn}">−</button>
        <button id="z_res" style="{btn};padding:0 10px;width:auto;">Reset</button>
      </div>
      <div id="treebox" style="width:100%;height:100%;cursor:grab;overflow:hidden;">{svg_text}</div>
    </div>
    <script>
    (function() {{
      var box = document.getElementById('treebox');
      var svg = box.querySelector('svg');
      if (!svg) return;
      svg.style.transformOrigin = '0 0';
      svg.style.width = '100%';
      svg.style.height = '100%';
      var scale = 1, tx = 0, ty = 0;
      function apply() {{
        svg.style.transform = 'translate(' + tx + 'px,' + ty + 'px) scale(' + scale + ')';
      }}
      function zoomAt(mx, my, factor) {{
        var ns = Math.min(40, Math.max(0.3, scale * factor));
        var k = ns / scale;
        tx = mx - k * (mx - tx);
        ty = my - k * (my - ty);
        scale = ns;
        apply();
      }}
      box.addEventListener('wheel', function(e) {{
        e.preventDefault();
        var r = box.getBoundingClientRect();
        zoomAt(e.clientX - r.left, e.clientY - r.top, e.deltaY < 0 ? 1.15 : 1 / 1.15);
      }}, {{ passive: false }});
      var drag = false, sx = 0, sy = 0;
      box.addEventListener('mousedown', function(e) {{
        drag = true; sx = e.clientX - tx; sy = e.clientY - ty; box.style.cursor = 'grabbing';
      }});
      window.addEventListener('mousemove', function(e) {{
        if (!drag) return; tx = e.clientX - sx; ty = e.clientY - sy; apply();
      }});
      window.addEventListener('mouseup', function() {{ drag = false; box.style.cursor = 'grab'; }});
      function centerZoom(factor) {{
        var r = box.getBoundingClientRect();
        zoomAt(r.width / 2, r.height / 2, factor);
      }}
      document.getElementById('z_in').onclick  = function() {{ centerZoom(1.25); }};
      document.getElementById('z_out').onclick = function() {{ centerZoom(1 / 1.25); }};
      document.getElementById('z_res').onclick = function() {{ scale = 1; tx = 0; ty = 0; apply(); }};
      apply();
    }})();
    </script>
    """
    components.html(html, height=height + 20)


def main():
    st.markdown(CSS, unsafe_allow_html=True)
    st.title("Baza homologów zuzalizyny")
    st.caption("Zuzalizyna (Q7MTD8, *Porphyromonas gingivalis*) — przegląd genomów "
               "bakterii i ich paralogów: domeny DUF, motyw katalityczny "
               "HExxHxxGxxH, peptyd sygnałowy i kontekst systemu sekrecji T9SS.")

    with st.sidebar:
        st.header("Dane")
        up = st.file_uploader("Wgraj swój plik CSV", type=["csv"])
        if up is not None:
            st.success("Używam wgranego pliku CSV.")

    if up is None:
        st.info("Wgraj plik CSV w panelu bocznym, aby zobaczyć bazę.")
        st.stop()

    raw = up.getvalue()

    try:
        df, genomes = process(raw)
    except Exception as e:
        st.error(f"Nie udało się przetworzyć pliku: {e}")
        st.stop()

    with st.sidebar:
        st.header("Filtry")
        genera = sorted(genomes["Genus"].dropna().unique().tolist())
        sel_genus = st.multiselect("Rodzaj", genera, default=[])

        t9_opts = [s for s in T9SS_ORDER if s in genomes.get("T9SS_status", pd.Series()).unique()]
        t9_opts += [s for s in genomes.get("T9SS_status", pd.Series()).dropna().unique()
                    if s not in t9_opts]
        sel_t9 = st.multiselect("Status T9SS", t9_opts, default=[])

        pmin, pmax = int(genomes["n_paralogs"].min()), int(genomes["n_paralogs"].max())
        if pmin < pmax:
            rng_par = st.slider("Liczba paralogów (genom)", pmin, pmax, (pmin, pmax))
        else:
            rng_par = (pmin, pmax)

        if "Per. ident" in genomes and genomes["Per. ident"].notna().any():
            imin = float(np.floor(genomes["Per. ident"].min()))
            imax = float(np.ceil(genomes["Per. ident"].max()))
            rng_id = st.slider("% identyczności do zapytania", imin, imax, (imin, imax))
        else:
            rng_id = None

        st.markdown("**Filtry paralogowe**")
        f_motif = st.selectbox("Motyw HExxHxxGxxH", ["wszystkie", "tylko z motywem",
                                                     "tylko bez motywu"])
        dom_filter = st.multiselect("Musi mieć domeny",
                                    ["DUF4953", "DUF5117", "DUF5118"], default=[])
        pred_opts = sorted(df["Prediction"].dropna().unique().tolist()) \
            if "Prediction" in df.columns else []
        sel_pred = st.multiselect("Predykcja peptydu sygn. (SignalP)", pred_opts, default=[])

    gmask = pd.Series(True, index=genomes.index)
    if sel_genus:
        gmask &= genomes["Genus"].isin(sel_genus)
    if sel_t9 and "T9SS_status" in genomes:
        gmask &= genomes["T9SS_status"].isin(sel_t9)
    gmask &= genomes["n_paralogs"].between(rng_par[0], rng_par[1])
    if rng_id is not None:
        gmask &= genomes["Per. ident"].between(rng_id[0], rng_id[1]) | genomes["Per. ident"].isna()
    fgen = genomes[gmask].copy()
    allowed = set(fgen["__genome"])

    pmask = df["__genome"].isin(allowed)
    if f_motif == "tylko z motywem":
        pmask &= df["has_HEXXH"]
    elif f_motif == "tylko bez motywu":
        pmask &= ~df["has_HEXXH"]
    for d in dom_filter:
        pmask &= df[f"has_{d}"]
    if sel_pred:
        pmask &= df["Prediction"].isin(sel_pred)
    fdf = df[pmask].copy()

    with st.sidebar:
        st.markdown("---")
        st.markdown(f"<span class='small'>Wynik filtrów: <b>{len(fgen)}</b> genomów, "
                    f"<b>{len(fdf)}</b> paralogów</span>", unsafe_allow_html=True)

    if fgen.empty:
        st.warning("Żaden genom nie spełnia wybranych filtrów.")
        st.stop()

    tab_dash, tab_browse, tab_par, tab_seq, tab_drzewo, tab_info = st.tabs(
        ["Statystyki", "Przeglądaj genomy", "Eksplorator paralogów",
         "Sekwencje", "Drzewo filogenetyczne", "Model danych"])

    with tab_dash:
        n_gen = len(fgen)
        n_par = len(fdf)
        mean_par = fgen["n_paralogs"].mean()
        top_row = fgen.loc[fgen["n_paralogs"].idxmax()]
        pct_complete = 100 * (fgen.get("T9SS_status") == "complete").mean() \
            if "T9SS_status" in fgen else np.nan
        pct_motif = 100 * fdf["has_HEXXH"].mean() if n_par else 0
        pct_full = 100 * fdf["has_all_DUF"].mean() if n_par else 0
        n_candidates = int(((fdf["architecture"] == ARCH_ORDER[0])
                            & (fdf["__genome"].map(
                                fgen.set_index("__genome").get("T9SS_status", pd.Series())
                                .to_dict()) == "complete")).sum()) \
            if "T9SS_status" in fgen else 0

        c = st.columns(4, gap="medium")
        kpi(c[0], "Genomy", fmt(n_gen, 0))
        kpi(c[1], "Paralogi", fmt(n_par, 0), f"średnio {mean_par:.1f} / genom")
        kpi(c[2], "Najwięcej paralogów", fmt(top_row["n_paralogs"], 0),
            str(top_row.get("Scientific Name", "")))
        kpi(c[3], "Genomy z kompletnym T9SS",
            "—" if pd.isna(pct_complete) else f"{pct_complete:.0f}%")
        c = st.columns(4, gap="medium")
        kpi(c[0], "Paralogi z motywem HExxHxxGxxH", f"{pct_motif:.0f}%")
        kpi(c[1], "Paralogi „pełne domenowo”", f"{pct_full:.0f}%",
            "wszystkie 3 domeny DUF")
        kpi(c[2], "Kompletny T9SS i SP/LIPO", fmt(n_candidates, 0),
            "peptyd + T9SS")
        kpi(c[3], "Rodzaje (genera)", fmt(fgen["Genus"].nunique(), 0))

        st.markdown("### Przegląd")
        a, b = st.columns(2)
        with a:
            vc = fgen["n_paralogs"].value_counts().sort_index().reset_index()
            vc.columns = ["Liczba paralogów", "Liczba genomów"]
            fig = px.bar(vc, x="Liczba paralogów", y="Liczba genomów",
                         template=PLOTLY_TEMPLATE, color_discrete_sequence=[PALETTE[0]],
                         title="Rozkład liczby paralogów na genom")
            fig.update_layout(bargap=0.12, height=360, margin=dict(t=50, b=10))
            st.plotly_chart(fig, use_container_width=True)
        with b:
            top = fgen.nlargest(15, "n_paralogs")
            fig = px.bar(top, x="n_paralogs", y="genome_label", orientation="h",
                         template=PLOTLY_TEMPLATE,
                         color="Genus", color_discrete_sequence=PALETTE,
                         title="Top 15 genomów wg liczby paralogów",
                         labels={"n_paralogs": "Paralogi", "genome_label": ""})
            fig.update_layout(height=360, margin=dict(t=50, b=10),
                              yaxis={"categoryorder": "total ascending"})
            st.plotly_chart(fig, use_container_width=True)

        tm = fgen.groupby(["Genus", "Scientific Name"], dropna=False)["n_paralogs"]\
            .sum().reset_index()
        fig = px.treemap(tm, path=[px.Constant("Wszystkie"), "Genus", "Scientific Name"],
                         values="n_paralogs", template=PLOTLY_TEMPLATE,
                         color="n_paralogs", color_continuous_scale="Teal",
                         title="Mapa drzewa: rodzaj → gatunek (rozmiar = liczba paralogów)")
        fig.update_layout(height=420, margin=dict(t=50, b=10))
        st.plotly_chart(fig, use_container_width=True)

        st.markdown("### System sekrecji T9SS")
        a, b = st.columns([1, 1])
        with a:
            if "T9SS_status" in fgen:
                sc = fgen["T9SS_status"].value_counts().reset_index()
                sc.columns = ["status", "n"]
                fig = px.pie(sc, names="status", values="n", hole=0.5,
                             template=PLOTLY_TEMPLATE,
                             color="status", color_discrete_map=T9SS_COLORS,
                             title="Status T9SS (genomy)")
                fig.update_layout(height=360, margin=dict(t=50, b=10))
                st.plotly_chart(fig, use_container_width=True)
        with b:
            comp_rows = []
            for comp in ["SprA", "PorN", "PorU", "PorV"]:
                col = f"{comp}_present"
                if col in fgen.columns:
                    comp_rows.append({"komponent": comp,
                                      "% genomów": 100 * (fgen[col] == "TAK").mean()})
            if comp_rows:
                cdf = pd.DataFrame(comp_rows)
                fig = px.bar(cdf, x="komponent", y="% genomów", template=PLOTLY_TEMPLATE,
                             color_discrete_sequence=[PALETTE[2]],
                             title="Obecność komponentów T9SS (SprA, PorN, PorU, PorV)")
                fig.update_layout(height=360, margin=dict(t=50, b=10), yaxis_range=[0, 100])
                st.plotly_chart(fig, use_container_width=True)

        st.markdown("### Domeny i motyw katalityczny (poziom paralogu)")
        a, b = st.columns(2)
        with a:
            dom_rows = []
            for name, col in [("DUF4953", "has_DUF4953"), ("DUF5117", "has_DUF5117"),
                              ("DUF5118", "has_DUF5118"), ("motyw HExxHxxGxxH", "has_HEXXH")]:
                dom_rows.append({"cecha": name, "% paralogów": 100 * fdf[col].mean()
                                 if len(fdf) else 0})
            ddf = pd.DataFrame(dom_rows)
            fig = px.bar(ddf, x="% paralogów", y="cecha", orientation="h",
                         template=PLOTLY_TEMPLATE, color_discrete_sequence=[PALETTE[0]],
                         title="Obecność domen i motywu")
            fig.update_layout(height=320, margin=dict(t=50, b=10), xaxis_range=[0, 100],
                              yaxis={"categoryorder": "total ascending"})
            st.plotly_chart(fig, use_container_width=True)
        with b:
            ac = fdf["architecture"].value_counts().reindex(ARCH_ORDER).dropna().reset_index()
            ac.columns = ["architecture", "n"]
            fig = px.pie(ac, names="architecture", values="n", hole=0.45,
                         template=PLOTLY_TEMPLATE, color="architecture",
                         color_discrete_map=ARCH_COLORS,
                         title="Klasy budowy paralogów")
            fig.update_layout(height=320, margin=dict(t=50, b=10),
                              legend=dict(orientation="h", y=-0.2))
            st.plotly_chart(fig, use_container_width=True)

        score_cols = [c for c in ["DUF4953_score", "DUF5117_score", "DUF5118_score"]
                      if c in fdf.columns]
        if score_cols:
            melt = fdf[score_cols].melt(var_name="domena", value_name="score").dropna()
            melt["domena"] = melt["domena"].str.replace("_score", "", regex=False)
            fig = px.box(melt, x="domena", y="score", color="domena", points="outliers",
                         template=PLOTLY_TEMPLATE, color_discrete_sequence=PALETTE,
                         title="Rozkład wyników (score) dopasowania domen HMM")
            fig.update_layout(height=360, margin=dict(t=50, b=10), showlegend=False)
            st.plotly_chart(fig, use_container_width=True)

        st.markdown("### Jakość i właściwości genomów")
        a, b = st.columns(2)
        with a:
            if "Per. ident" in fgen:
                fig = px.histogram(fgen, x="Per. ident", nbins=25, template=PLOTLY_TEMPLATE,
                                   color_discrete_sequence=[PALETTE[1]],
                                   title="Rozkład % identyczności do zuzalizyny")
                fig.update_layout(height=340, margin=dict(t=50, b=10))
                st.plotly_chart(fig, use_container_width=True)
        with b:
            if {"Assembly Stats GC Percent"}.issubset(fgen.columns):
                fig = px.scatter(fgen, x="Assembly Stats GC Percent", y="n_paralogs",
                                 color="Genus", size="n_paralogs", template=PLOTLY_TEMPLATE,
                                 color_discrete_sequence=PALETTE,
                                 hover_name="genome_label",
                                 title="GC% genomu a liczba paralogów",
                                 labels={"Assembly Stats GC Percent": "GC %",
                                         "n_paralogs": "Paralogi"})
                fig.update_layout(height=340, margin=dict(t=50, b=10))
                st.plotly_chart(fig, use_container_width=True)

        if {"CheckM completeness", "CheckM contamination"}.issubset(fgen.columns):
            fig = px.scatter(fgen, x="CheckM completeness", y="CheckM contamination",
                             size="n_paralogs", color="Genus", template=PLOTLY_TEMPLATE,
                             color_discrete_sequence=PALETTE, hover_name="genome_label",
                             title="Jakość genomu wg CheckM (kompletność vs kontaminacja)")
            fig.update_layout(height=380, margin=dict(t=50, b=10))
            st.plotly_chart(fig, use_container_width=True)

    with tab_browse:
        st.markdown("#### Lista genomów")
        sort_opts = {
            "Nazwa (A→Z)": ("Scientific Name", True),
            "Nazwa (Z→A)": ("Scientific Name", False),
            "Liczba paralogów ↓": ("n_paralogs", False),
            "Liczba paralogów ↑": ("n_paralogs", True),
            "% identyczności ↓": ("Per. ident", False),
            "GC % ↓": ("Assembly Stats GC Percent", False),
            "Kompletność CheckM ↓": ("CheckM completeness", False),
            "Paralogi z motywem ↓": ("n_HEXXH", False),
        }
        cc = st.columns([2, 3])
        sort_choice = cc[0].selectbox("Sortuj wg", list(sort_opts.keys()))
        sort_col, asc = sort_opts[sort_choice]
        if sort_col in fgen.columns:
            fgen_sorted = fgen.sort_values(sort_col, ascending=asc,
                                           na_position="last").reset_index(drop=True)
        else:
            fgen_sorted = fgen.reset_index(drop=True)

        show_cols = ["genome_label", "Genus", "n_paralogs", "T9SS_status", "Per. ident",
                     "Assembly Stats GC Percent", "CheckM completeness", "n_HEXXH",
                     "n_full", "n_secreted", "Assembly Level", "Assembly"]
        show_cols = [c for c in show_cols if c in fgen_sorted.columns]
        st.dataframe(
            fgen_sorted[show_cols], use_container_width=True, hide_index=True, height=330,
            column_config={
                "genome_label": st.column_config.TextColumn("Organizm", width="large"),
                "Genus": "Rodzaj",
                "n_paralogs": st.column_config.NumberColumn("Paralogi", format="%d"),
                "T9SS_status": "T9SS",
                "Per. ident": st.column_config.ProgressColumn(
                    "% ident.", min_value=0, max_value=100, format="%.1f"),
                "Assembly Stats GC Percent": st.column_config.NumberColumn("GC %", format="%.1f"),
                "CheckM completeness": st.column_config.ProgressColumn(
                    "CheckM %", min_value=0, max_value=100, format="%.0f"),
                "n_HEXXH": st.column_config.NumberColumn("z motywem", format="%d"),
                "n_full": st.column_config.NumberColumn("pełne dom.", format="%d"),
                "n_secreted": st.column_config.NumberColumn("z peptydem", format="%d"),
                "Assembly Level": "Poziom",
            })

        labels = fgen_sorted["genome_label"] + "  ·  " + \
            fgen_sorted["n_paralogs"].astype(int).astype(str) + " paral."
        pick = st.selectbox("Wybierz genom do szczegółów", labels.tolist())
        gidx = labels.tolist().index(pick)
        grow = fgen_sorted.iloc[gidx]
        gkey = grow["__genome"]
        paras = df[df["__genome"] == gkey].copy().sort_values("paralog_num")

        st.markdown(f"### {grow.get('Scientific Name','')}")
        meta = st.columns([3, 2])
        with meta[0]:
            url = grow.get("Accession_URL")
            link = f" · [NCBI ↗]({url})" if isinstance(url, str) else ""
            st.markdown(
                f"<div class='detailcard'>"
                f"<b>Szczep:</b> {grow.get('Organism Infraspecific Names Strain','—')} &nbsp;|&nbsp; "
                f"<b>Taxid:</b> {fmt(grow.get('Taxid'),0)}<br>"
                f"<b>Assembly:</b> {grow.get('Assembly','—')} ({grow.get('Assembly Name','—')})<br>"
                f"<b>Accession:</b> {grow.get('Accession_ID','—')}{link}<br>"
                f"<b>BioProject / BioSample:</b> {grow.get('Assembly BioProject Accession','—')} / "
                f"{grow.get('Assembly BioSample Accession','—')}<br>"
                f"<b>Sekwencjonowanie:</b> {grow.get('Assembly Sequencing Tech','—')} &nbsp;|&nbsp; "
                f"<b>Materiał typowy:</b> {grow.get('Type Material Display Text','—')}"
                f"</div>", unsafe_allow_html=True)
        with meta[1]:
            t9 = grow.get("T9SS_status", "—")
            color = T9SS_COLORS.get(t9, "#6b7280")
            st.markdown(f"<div class='detailcard'><b>System sekrecji T9SS</b><br>"
                        f"<span class='badge' style='background:{color}'>{t9}</span> "
                        f"&nbsp; znaleziono komponentów: "
                        f"<b>{fmt(grow.get('T9SS_components_found'),0)}</b></div>",
                        unsafe_allow_html=True)

        m = st.columns(5, gap="medium")
        kpi(m[0], "Paralogi", fmt(grow["n_paralogs"], 0))
        kpi(m[1], "Wielkość genomu", f"{fmt(grow.get('Assembly Stats Total Sequence Length',np.nan)/1e6,2)} Mb"
            if pd.notna(grow.get('Assembly Stats Total Sequence Length', np.nan)) else "—")
        kpi(m[2], "GC %", fmt(grow.get("Assembly Stats GC Percent"), 1))
        kpi(m[3], "CheckM kompletność", f"{fmt(grow.get('CheckM completeness'),1)}%")
        kpi(m[4], "CheckM kontaminacja", f"{fmt(grow.get('CheckM contamination'),1)}%")

        comp_tab = []
        for comp in ["SprA", "PorN", "PorU", "PorV"]:
            if f"{comp}_present" in grow.index:
                comp_tab.append({"Komponent": comp,
                                 "Obecny": yn(grow[f"{comp}_present"] == "TAK"),
                                 "e-value": grow.get(f"{comp}_evalue"),
                                 "Liczba kopii": grow.get(f"{comp}_count")})
        if comp_tab:
            st.markdown("**Komponenty T9SS:**")
            st.dataframe(pd.DataFrame(comp_tab), hide_index=True, use_container_width=True)

        st.markdown(f"#### Paralogi tego genomu ({len(paras)})")
        pcols = ["Paralog_ID", "ID_bialka", "paralog_num", "length_aa", "Prediction",
                 "CS Position", "czy_jest_HEXXHXXGXXH", "czy_jest_DUF4953",
                 "czy_jest_DUF5117", "czy_jest_DUF5118", "DUF4953_score", "DUF5117_score",
                 "DUF5118_score", "completeness_score", "architecture"]
        pcols = [c for c in pcols if c in paras.columns]
        ptab = styled_bool_df(paras[pcols], BOOL_COLS)
        st.dataframe(
            ptab, hide_index=True, use_container_width=True,
            column_config={
                "Paralog_ID": st.column_config.TextColumn("ID paralogu", width="medium"),
                "ID_bialka": st.column_config.TextColumn("ID białka (NCBI)", width="medium"),
                "paralog_num": st.column_config.NumberColumn("#", format="%d"),
                "length_aa": st.column_config.NumberColumn("dł. (aa)", format="%d"),
                "Prediction": "Peptyd sygn.",
                "czy_jest_HEXXHXXGXXH": "motyw",
                "czy_jest_DUF4953": "DUF4953", "czy_jest_DUF5117": "DUF5117",
                "czy_jest_DUF5118": "DUF5118",
                "DUF4953_score": st.column_config.NumberColumn("4953 sc.", format="%.0f"),
                "DUF5117_score": st.column_config.NumberColumn("5117 sc.", format="%.0f"),
                "DUF5118_score": st.column_config.NumberColumn("5118 sc.", format="%.0f"),
                "completeness_score": st.column_config.ProgressColumn(
                    "kompletność /5", min_value=0, max_value=5, format="%d"),
                "architecture": st.column_config.TextColumn("budowa", width="medium"),
            })
        st.caption("Długość (aa) liczona z pełnej sekwencji, gdy jest dostępna "
                   "(kolumna `sekwencje_cale`); w przeciwnym razie z dopasowania BLAST.")
        st.download_button("Pobierz paralogi tego genomu (FASTA)",
                           to_fasta(paras), file_name=f"{gkey}_paralogi.fasta")

    with tab_par:
        st.markdown("#### Wszystkie paralogi (po filtrach)")
        q = st.text_input("Szukaj (gatunek / Paralog_ID / ID białka / Assembly)", "")
        view = fdf.copy()
        if q:
            ql = q.lower()
            id_bialka = (view["ID_bialka"] if "ID_bialka" in view.columns
                         else pd.Series("", index=view.index))
            hay = (view["Scientific Name"].fillna("") + " "
                   + view["Paralog_ID"].fillna("") + " "
                   + id_bialka.fillna("") + " "
                   + view["__genome"].fillna("")).str.lower()
            view = view[hay.str.contains(ql, regex=False)]

        cols = ["Scientific Name", "Genus", "Paralog_ID", "ID_bialka", "length_aa",
                "Per. ident", "T9SS_status", "Prediction", "CS Position",
                "czy_jest_HEXXHXXGXXH", "czy_jest_DUF4953", "czy_jest_DUF5117",
                "czy_jest_DUF5118", "completeness_score", "architecture", "Assembly"]
        cols = [c for c in cols if c in view.columns]
        st.caption(f"Wyświetlono {len(view)} paralogów.")
        st.dataframe(styled_bool_df(view[cols], BOOL_COLS), hide_index=True,
                     use_container_width=True, height=520,
                     column_config={
                         "Scientific Name": st.column_config.TextColumn("Organizm", width="medium"),
                         "ID_bialka": st.column_config.TextColumn("ID białka (NCBI)", width="medium"),
                         "length_aa": st.column_config.NumberColumn("dł. (aa)", format="%d"),
                         "Per. ident": st.column_config.NumberColumn("% ident.", format="%.1f"),
                         "completeness_score": st.column_config.ProgressColumn(
                             "kompl./5", min_value=0, max_value=5, format="%d"),
                     })
        d1, d2 = st.columns(2)
        d1.download_button("Pobierz tabelę paralogów (CSV)",
                           view[cols].to_csv(index=False).encode("utf-8"),
                           file_name="paralogi_filtr.csv")
        d2.download_button("Pobierz sekwencje (FASTA)", to_fasta(view),
                           file_name="paralogi_filtr.fasta")

    with tab_seq:
        st.markdown("#### Podgląd sekwencji z zaznaczeniem motywu i peptydu sygnałowego")
        g_labels = fgen.sort_values("Scientific Name")["genome_label"].tolist()
        gpick = st.selectbox("Genom", g_labels, key="seqgenome")
        gk = fgen[fgen["genome_label"] == gpick]["__genome"].iloc[0]
        sub = df[df["__genome"] == gk].sort_values("paralog_num")
        ppick = st.selectbox("Paralog", sub["Paralog_ID"].tolist())
        pr = sub[sub["Paralog_ID"] == ppick].iloc[0]

        # Wybór sekwencji: pełna ("sekwencje_cale") gdy jest, inaczej dopasowanie BLAST.
        full_seq = pr.get("sekwencje_cale")
        has_full = isinstance(full_seq, str) and full_seq.strip() != ""
        display_seq = full_seq if has_full else (pr.get("sequence") or "")

        cc = st.columns(4, gap="medium")
        kpi(cc[0], "Długość", f"{fmt(len(display_seq),0)} aa")
        kpi(cc[1], "Peptyd sygnałowy", pr.get("Prediction", "—"))
        kpi(cc[2], "Motyw HExxHxxGxxH", yn(pr["has_HEXXH"]))
        kpi(cc[3], "Domeny DUF (z 3)", fmt(pr["n_DUF"], 0))

        prot_id = pr.get("ID_bialka")
        if isinstance(prot_id, str) and prot_id.strip():
            st.markdown(f"<span class='small'>ID białka (NCBI): <b>{prot_id}</b></span>",
                        unsafe_allow_html=True)

        if not has_full:
            st.warning("UWAGA: to nie cała sekwencja, tylko dopasowanie")
            st.caption("Brak pełnej sekwencji (kolumna `sekwencje_cale`) — "
                       "peptyd sygnałowy nie jest zaznaczany, bo jego pozycja odnosi się "
                       "do pełnego białka, a nie do fragmentu dopasowania.")
            legend = ("<span class='small'>Legenda: "
                      "<span class='badge' style='background:#c2503a'>motyw katalityczny "
                      "HExxHxxGxxH</span></span>")
        else:
            legend = ("<span class='small'>Legenda: "
                      "<span class='badge' style='background:#345'>peptyd sygnałowy</span> "
                      "<span class='badge' style='background:#c2503a'>motyw katalityczny "
                      "HExxHxxGxxH</span></span>")
        st.markdown(legend, unsafe_allow_html=True)

        # peptyd sygnałowy zaznaczamy tylko, gdy mamy pełną sekwencję
        cs_first = parse_cs_first(pr.get("CS Position")) if has_full else None
        st.markdown(highlight_sequence_html(display_seq, cs_first),
                    unsafe_allow_html=True)
        motifs = [f"{m.start()+1}–{m.end()}" for m in MOTIF_REGEX.finditer(display_seq)]
        if motifs:
            st.caption("Pozycje motywu: " + ", ".join(motifs))
        st.download_button("Pobierz tę sekwencję (FASTA)",
                           to_fasta(sub[sub["Paralog_ID"] == ppick]),
                           file_name=f"{ppick}.fasta")

    with tab_drzewo:
        st.markdown("#### Drzewo filogenetyczne")
        st.caption("Wgraj plik SVG z drzewem — możesz je przybliżać kółkiem myszy, "
                   "przesuwać przeciąganiem i korzystać z przycisków +/−/reset.")
        tree_svg = st.file_uploader("Plik SVG z drzewem", type=["svg"], key="tree_svg")
        if tree_svg is not None:
            svg_text = tree_svg.getvalue().decode("utf-8", errors="replace")
            pokaz_drzewo_svg(svg_text, height=720)
            st.download_button("Pobierz SVG z drzewem", tree_svg.getvalue(),
                               file_name="drzewo_filogenetyczne.svg",
                               mime="image/svg+xml", key="tree_dl")
        else:
            st.info("Wgraj SVG z drzewem, aby je zobaczyć.")

    with tab_info:
        st.markdown("""


#### Sekwencje — co jest wyświetlane

* **`sekwencje_cale`** to **pełna sekwencja białka** (z N-końcem i peptydem
  sygnałowym). To ją aplikacja pokazuje domyślnie i na niej zaznacza **peptyd
  sygnałowy** (na podstawie `CS Position` z SignalP) oraz **motyw HExxHxxGxxH**.
* Jeżeli `sekwencje_cale` jest pusta, pokazywana jest **`sequence`** — czyli sam
  **fragment dopasowania BLAST**. Wtedy wyświetla się komunikat
  *„UWAGA: to nie cała sekwencja, tylko dopasowanie”*, a **peptyd sygnałowy nie jest
  zaznaczany** 
* **`ID_bialka`** to numer dostępu białka w NCBI (np. `WP_012458658.1`) — widoczny
  w tabelach paralogów, w podglądzie sekwencji i w nagłówkach FASTA.


#### Wskazówki
* Trzeba wgrać CSV z danymi w panelu bocznym — kolumny rozpoznawane są po nazwach,
  brakujące są pomijane bez błędu.
        """)
        with st.expander("Podgląd surowych danych (po przetworzeniu)"):
            st.dataframe(df.head(50), use_container_width=True)


if __name__ == "__main__":
    main()