"""o5c lane B: src-literature-138 (figshare 31083118, Starkweather et al., 'Rapid Homotopic Communication Between Human
Orbitofrontal Subregions') -> iEEG-BIDS DERIVATIVE (authors' preprocessed timeseries). Per subject folder (figshare file
ids, grouped as listed by the figshare API): condensed_ecog_OFC_only.mat (MAT v7.3 struct array, one element per OFC
electrode: unreferenced_ecog and bipolar timeseries at Fs = 1000 Hz, brain_region, contact, channel_names = reference,
sulcus-relative coordinates lat/med/mid/olf/trans_coor), stimset.mat (single-pulse stimulation times/amplitudes per
stimulating electrode), metadata.mat (CCEP measures per stim->rec pair; sourcedata), condensed_ecog_info_only.mat.
Two acquisitions per subject: acq-unreferenced and acq-bipolar (BrainVision float32, uV). Stimulation pulses -> events.tsv.
Usage: python o5c_B_conv_lit138.py src-literature-138"""
import glob, json, os, re, shutil, sys, urllib.request
import h5py, numpy as np, scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o5c_B_common import *  # noqa
from b3w10_common import write_run  # noqa

CID = sys.argv[1] if len(sys.argv) > 1 else "src-literature-138"
SRC = src_dir(CID); OUT = fresh_out(CID); REP = rep_dir(CID)
F = os.path.join(SRC, "files")
FIG = "31083118"; DOI = f"10.6084/m9.figshare.{FIG}.v2"; URL = f"https://figshare.com/articles/dataset/{FIG}"
CONV = "harness/o5c_B_conv_lit138.py (ieeg-nemar-20261005 campaign, lane o5c-B)"
# figshare file ids per subject folder, in the order the figshare API lists them (folder 'subject N' of the release scripts)
SUBJ = {"01": {"info": "61081318", "ofc": "61081396", "meta": "61081321", "stim": "61081324"},
        "02": {"info": "61081327", "ofc": "61081402", "meta": "61081330", "stim": "61081333"},
        "03": {"info": "61081336", "ofc": "61081393", "meta": "61081339", "stim": "61081345"},
        "04": {"info": "61081342", "ofc": "61081411", "meta": "61081348", "stim": "61081351"},
        "05": {"info": "61081354", "ofc": "61081408", "meta": "61081357", "stim": "61081360"},
        "06": {"info": "61081363", "ofc": "61081405", "meta": "61081366", "stim": "61081369"},
        "07": {"info": "61081372", "ofc": "61081399", "meta": "61081375", "stim": "61081378"}}
try:
    rec = json.load(urllib.request.urlopen(f"https://api.figshare.com/v2/articles/{FIG}", timeout=60))
except Exception as e:
    rec = {"error": repr(e)}
SD = os.path.join(OUT, "sourcedata", f"figshare-{FIG}"); os.makedirs(SD, exist_ok=True)
wjson(os.path.join(SD, f"figshare_article_{FIG}.json"), rec)
authors = [a["full_name"] for a in rec.get("authors", [])] or ["Clara Starkweather"]
lic = (rec.get("license") or {}).get("name", "CC BY 4.0")
rep = {"cand_id": CID, "converter": CONV, "licence": {"value": lic, "read_from": "figshare api license.name"}, "runs": [], "sourcedata": [], "subjects": {}}
COORD = ["lat_coor", "med_coor", "mid_coor", "olf_coor", "trans_coor"]


def h5str(f, r):
    o = f[r]
    return "".join(chr(int(c)) for c in np.asarray(o[()]).flatten(order="F"))


def load_struct(p, name):
    """MAT v7.3 struct array -> dict field -> list of per-element values (strings, scalars or 1-D arrays)."""
    out = {}
    with h5py.File(p, "r") as f:
        g = f[name]
        for k in g.keys():
            ds = g[k]
            vals = []
            if ds.dtype.kind == "O":
                for r in np.asarray(ds[()]).flatten():
                    o = f[r]
                    mc = o.attrs.get("MATLAB_class", b"")
                    mc = mc.decode() if isinstance(mc, bytes) else str(mc)
                    if mc == "char":
                        vals.append(h5str(f, r))
                    else:
                        a = np.asarray(o[()])
                        vals.append(float(a.flatten()[0]) if a.size == 1 else a.flatten().astype(np.float64))
            else:
                a = np.asarray(ds[()]).flatten()
                vals = [float(x) for x in a]
            out[k] = vals
    return out


scans = {}
for sub, ids in SUBJ.items():
    ofc = os.path.join(F, ids["ofc"], "condensed_ecog_OFC_only.mat")
    S = load_struct(ofc, "condensed_ecog")
    n = len(S["brain_region"])
    fs = float(S["Fs"][0]) if "Fs" in S else 1000.0
    fs_src = "Fs field" if "Fs" in S else "assumed 1000 Hz (Fs field absent in this file; all other files state 1000; release scripts use 1000)"
    regions = S["brain_region"]; contacts = [int(c) for c in S["contact"]]
    names = [f"{r}{c}" for r, c in zip(regions, contacts)]
    if len(set(names)) != n:  # disambiguate with the brain_region_number
        names = [f"{r}{int(S['brain_region_number'][i])}c{c}" for i, (r, c) in enumerate(zip(regions, contacts))]
    assert len(set(names)) == n, names
    ref = S.get("channel_names", ["n/a"])[0] if "channel_names" in S else "n/a"
    D = os.path.join(OUT, f"sub-{sub}", "ieeg")
    info = {"n_channels": n, "fs": fs, "fs_source": fs_src, "reference_channel": ref, "regions": sorted(set(regions)), "acq": {}}
    for field, acq in (("unreferenced_ecog", "unreferenced"), ("bipolar", "bipolar")):
        if field not in S:
            info["acq"][acq] = "field absent"; continue
        cols = S[field]
        lens = [len(np.atleast_1d(c)) for c in cols]
        L = max(lens)
        X = np.full((L, n), np.nan, dtype=np.float64)
        bad = []
        for i, c in enumerate(cols):
            c = np.atleast_1d(np.asarray(c, dtype=np.float64))
            if c.size < 2 or not np.isfinite(c).any():
                bad.append(names[i]); continue
            X[:c.size, i] = c
        fin = X[np.isfinite(X)]
        absmax = float(np.max(np.abs(fin))) if fin.size else 0.0
        scale_note = "values stored in volts in the release (|x| < 1); multiplied by 1e6" if absmax < 1 else "values stored in microvolts in the release"
        if absmax < 1:
            X = X * 1e6
        stem = f"sub-{sub}_task-ccep_acq-{acq}"
        ch_extra = [{"status": "bad" if nm in bad else "good", "status_description": ("no signal in the release (empty or all-NaN; first contact of a bipolar chain)" if nm in bad else "n/a"),
                     "group": re.sub(r"\d+$", "", regions[i]), "reference": ("bipolar (adjacent contact)" if acq == "bipolar" else (ref if ref != "n/a" else "n/a")),
                     "description": f"OFC electrode {regions[i]} contact {contacts[i]} (release struct element {i + 1}); authors' preprocessed {field} timeseries; {scale_note}"} for i, nm in enumerate(names)]
        # events: stimulation pulses
        ev = []
        stim = os.path.join(F, ids["stim"], "stimset.mat")
        if os.path.exists(stim):
            m = sio.loadmat(stim, squeeze_me=False, struct_as_record=True)
            el = m["stimset"][0, 0]["electrode"]
            for e in range(el.shape[1]):
                it = el[0, e]["intensity"]
                if not it.size:
                    continue
                for j in range(it.shape[1]):
                    amp = np.asarray(it[0, j]["amplitude"]).flatten(); tm = np.asarray(it[0, j]["times"]).flatten()
                    a = float(amp[0]) if amp.size else "n/a"
                    for t in tm:
                        ev.append([round(float(t) / fs, 6), "n/a", int(t), f"stim_{names[e] if e < n else 'electrode' + str(e + 1)}", a, e + 1])
            ev.sort(key=lambda r: r[0])
        r = write_run(D, stem, X, fs, names, ["µV"] * n, ["ECOG"] * n, {
            "TaskName": "ccep", "TaskDescription": "Cortico-cortical evoked potential (CCEP) mapping: single-pulse direct electrical stimulation of orbitofrontal electrodes while recording from the other orbitofrontal contacts (release description).",
            "PowerLineFrequency": 60, "SoftwareFilters": {"authors' preprocessing": {"description": "condensed/preprocessed timeseries as released (details in the article); bipolar = adjacent-contact re-referencing by the authors"}},
            "HardwareFilters": "n/a", "iEEGReference": "bipolar (adjacent contacts; authors)" if acq == "bipolar" else (f"recording reference channel {ref} (release field channel_names)" if ref != "n/a" else "n/a"),
            "iEEGPlacementScheme": f"Orbitofrontal cortex electrodes (regions in the release: {', '.join(sorted(set(regions)))}); electrode types are not stated in the release (ECOG used for all)",
            "ElectricalStimulation": True, "ElectricalStimulationParameters": "Single-pulse stimulation; amplitudes (mA, release field amplitude) and pulse times per stimulating electrode in events.tsv", "Manufacturer": "n/a", "InstitutionName": "n/a",
            "RecordingType": "continuous"}, ch_extra=ch_extra, comment=f"{CONV}: {ids['ofc']}/condensed_ecog_OFC_only.mat field {field}")
        if ev:
            wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "sample", "trial_type", "stim_amplitude_mA", "stim_electrode_index"], ev)
            wjson(os.path.join(D, stem + "_events.json"), {"trial_type": {"Description": "stim_<channel> = single stimulation pulse delivered at that electrode (release stimset.mat)"},
                  "stim_amplitude_mA": {"Description": "stimulation amplitude of the pulse train (release field amplitude; unit assumed mA as in the article)", "Units": "mA"},
                  "stim_electrode_index": {"Description": "1-based index of the stimulating electrode in the release struct (same indexing as condensed_ecog)"},
                  "sample": {"Description": "release field times (sample index at 1000 Hz)"}})
        r.update(sub=sub, acq=acq, stem=stem, bad=bad, absmax=absmax, scale_note=scale_note, n_events=len(ev), lens=sorted(set(lens)))
        rep["runs"].append(r); info["acq"][acq] = {"n_bad": len(bad), "duration_s": r["duration_s"], "roundtrip_ok": r["roundtrip_ok"]}
        scans.setdefault(sub, []).append([f"ieeg/{stem}_ieeg.vhdr", "n/a"])
        print(stem, json.dumps({k: r[k] for k in ("n_samples", "n_channels", "roundtrip_ok", "n_events", "scale_note")}, default=str), flush=True)
        del X
    # electrodes: sulcus-relative coordinates as extra columns
    rows = []
    for i, nm in enumerate(names):
        rows.append([nm, "n/a", "n/a", "n/a", "n/a", re.sub(r"\d+$", "", regions[i]), regions[i], S.get("anatomical_info", ["n/a"] * n)[i]] + [S[c][i] if c in S else "n/a" for c in COORD])
    wtsv(os.path.join(D, f"sub-{sub}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "group", "brain_region", "anatomical_info"] + COORD, rows)
    wjson(os.path.join(D, f"sub-{sub}_space-Other_coordsystem.json"), {"iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "mm",
          "iEEGCoordinateSystemDescription": "No x/y/z coordinates are released. The release gives per-contact distances relative to orbitofrontal sulci (columns lat_coor, med_coor, mid_coor, olf_coor, trans_coor, in mm): the release description defines lat_coor as the medial/lateral distance from the lateral orbital sulcus; the other columns are copied with their release names (presumably distances from the medial orbital sulcus, midline, olfactory sulcus and transverse orbital sulcus; not defined in the release). brain_region = release electrode label."})
    rep["subjects"][sub] = info
    # sourcedata: the four release files of this subject, verbatim
    for k, fid in ids.items():
        for p in glob.glob(os.path.join(F, fid, "*")):
            dst = os.path.join(SD, f"subject_{sub}", os.path.basename(p))
            c = copy_verbatim(p, dst); rep["sourcedata"].append({"dst": os.path.relpath(dst, OUT), "figshare_file_id": fid, "bytes": c["bytes"], "sha256": c["sha256"]})
for p in glob.glob(os.path.join(F, "*.m")):
    c = copy_verbatim(p, os.path.join(SD, "code", os.path.basename(p))); rep["sourcedata"].append({"dst": os.path.relpath(c["dst"], OUT), "bytes": c["bytes"], "sha256": c["sha256"]})
for sub, rows in scans.items():
    wtsv(os.path.join(OUT, f"sub-{sub}", f"sub-{sub}_scans.tsv"), ["filename", "acq_time"], rows)
wjson(os.path.join(OUT, "space-Other_electrodes.json"), {"SpatialReference": "n/a: no imaging is released; the columns lat/med/mid/olf/trans_coor are distances (mm) relative to orbitofrontal sulci as given by the authors, not image coordinates", "brain_region": {"Description": "release electrode label (e.g. ROFC = right orbitofrontal, LOFC = left orbitofrontal, OL/OD/LOD/ROD/RSFG/RSGC/LSFG as named by the authors; abbreviations not expanded)"},
      "anatomical_info": {"Description": "release field anatomical_info ('orb' = orbital)"}, "lat_coor": {"Description": "medial/lateral distance from the lateral orbital sulcus (release description)", "Units": "mm"},
      "med_coor": {"Description": "release field med_coor (not defined in the release)", "Units": "mm"}, "mid_coor": {"Description": "release field mid_coor (not defined in the release)", "Units": "mm"},
      "olf_coor": {"Description": "release field olf_coor (not defined in the release)", "Units": "mm"}, "trans_coor": {"Description": "release field trans_coor (not defined in the release)", "Units": "mm"}})
subs = sorted(SUBJ)
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "group", "release_folder"], [[f"sub-{s}", "n/a", "n/a", "epilepsy", f"subject {int(s)} (figshare files {', '.join(SUBJ[s].values())})"] for s in subs])
wjson(os.path.join(OUT, "participants.json"), {"age": {"Description": "not released"}, "sex": {"Description": "not released"}, "group": {"Description": "patients with intracranial electrodes for epilepsy monitoring (CCEP mapping study; not stated per participant in the release)"},
      "release_folder": {"Description": "release folder name used by the analysis scripts and the figshare file ids of its files"}})
hours = sum(r["duration_s"] for r in rep["runs"] if r["acq"] == "unreferenced") / 3600
wjson(os.path.join(OUT, "dataset_description.json"), dataset_description(
    "StarkweatherOFC-CCEP: orbitofrontal cortico-cortical evoked potential recordings, 7 patients (Starkweather et al.) (derivative)",
    authors, "CC-BY-4.0", DOI, URL, [f"https://doi.org/{DOI}"],
    f"Cite the figshare dataset doi:{DOI} and the manuscript \"Rapid Homotopic Communication Between Human Orbitofrontal Subregions\" (Starkweather et al.).",
    ["Collected under the source study's approval; the figshare record does not state the committee (see the manuscript)"],
    [{"Name": "o5c_B_conv_lit138.py", "Description": "Authors' preprocessed per-electrode timeseries (unreferenced and bipolar) repackaged as BrainVision float32 in microvolts; stimulation pulses from stimset.mat to events.tsv; all release files in sourcedata/. No further processing. Code in code/."}],
    dataset_type="derivative", keywords=["orbitofrontal cortex", "cortico-cortical evoked potentials", "CCEP", "single-pulse electrical stimulation", "intracranial EEG", "effective connectivity"],
    source_version=f"figshare version {rec.get('version', 'n/a')} ({rec.get('published_date', 'n/a')})"))
desc = re.sub("<[^>]+>", " ", rec.get("description", "")).strip()
lines = "\n".join(f"  - sub-{s}: {rep['subjects'][s]['n_channels']} OFC contacts ({', '.join(rep['subjects'][s]['regions'])}), {rep['subjects'][s]['acq'].get('unreferenced', {}).get('duration_s', 0) / 60:.1f} min, reference {rep['subjects'][s]['reference_channel']}" for s in subs)
open(os.path.join(OUT, "README"), "w").write(f"""# Orbitofrontal cortico-cortical evoked potential (CCEP) recordings, 7 patients (Starkweather et al.) - derivative

Release: figshare {FIG} "{rec.get('title', '')}" ({', '.join(authors)}; doi:{DOI}, {lic}). Release description:
"{desc}"

## Content (7 participants, {hours:.2f} h per acquisition)
For each participant the release holds `condensed_ecog_OFC_only.mat`, a struct array with one element per orbitofrontal
electrode contact carrying the authors' preprocessed timeseries at 1000 Hz: `unreferenced_ecog` (written here as
`acq-unreferenced`) and `bipolar` (adjacent-contact re-referenced, `acq-bipolar`; the first contact of a chain has no
bipolar signal and is marked `bad`). Values are written in microvolts (files whose release values are in volts are
multiplied by 1e6; stated per channel in channels.tsv). This is a DERIVATIVE dataset: the raw recordings are not in the
release. Channel names are `<brain_region><contact>` from the release fields. Single-pulse stimulation times and
amplitudes (`stimset.mat`) are in events.tsv (onset = release sample / 1000).
{lines}

## Localisation
No x/y/z coordinates are released. electrodes.tsv carries the release's sulcus-relative distances per contact
(`lat_coor` = medial/lateral distance from the lateral orbital sulcus per the release description; `med_coor`,
`mid_coor`, `olf_coor`, `trans_coor` copied with their release names) and the release electrode label (`brain_region`).

## sourcedata/figshare-{FIG}/
All release files per subject folder, byte for byte: `condensed_ecog_OFC_only.mat`, `condensed_ecog_info_only.mat` (the
files of subjects 04 and 07 are truncated in the release, 2.9 kB / 2.2 kB, and cannot be opened), `metadata.mat` (CCEP
onset latency / peak amplitude / peak latency / significance per stimulated -> recorded pair), `stimset.mat`, and the
authors' MATLAB scripts under `code/`.

## Provenance
Converted by the ieeg-nemar-20261005 campaign (lane o5c-B) with `code/o5c_B_conv_lit138.py`; round-trip checks (MNE
read-back) in code/conversion_report.json. A different dataset of the same first author (approach-avoidance SEEG HFA) is
NEMAR nm000383; participant overlap is unknown.
""")
copy_code(CID, ["o5c_B_conv_lit138.py", "o5c_B_common.py", "b3w10_common.py", "b3w3_common.py"])
bidsignore(OUT, ["sourcedata/", "code/"])
rep["hours"] = hours; rep["subjects_n"] = len(subs); rep["n_recordings"] = len(rep["runs"])
rep["roundtrip_all_ok"] = all(r["roundtrip_ok"] for r in rep["runs"])
report(CID, rep)
print(json.dumps({"runs": len(rep["runs"]), "hours": hours, "rt": rep["roundtrip_all_ok"]}, indent=1, default=str))
print("CONVERT_OK" if rep["roundtrip_all_ok"] else "CONVERT_ROUNDTRIP_FAIL")
