"""The plugin window: pick structures and parameters, run, inspect the result.

Three decisions shape this file.

**The analysis runs on a worker thread.** Building networks for a few hundred
structures takes minutes. On Qt's main thread that freezes PyMOL solid -- no
redraw, no rotation, and on Windows the title bar goes "(Not Responding)"
mid-demonstration. :class:`_Worker` does the computation;
:class:`WatConDialog` only draws.

**Only the main thread talks to PyMOL.** ``cmd`` is not safe to drive from a
worker. This is affordable precisely because :mod:`WatCon.scene` computes the
whole picture without importing PyMOL, so the thread boundary falls in a natural
place: the worker returns a ``Scene``, the main thread displays it.

**The session is displayed by running the same ``.pml`` that ``watcon view``
writes.** Not by replaying the command list separately -- by writing the file and
``@``-ing it. So the plugin and the command line cannot produce different
pictures; they execute the identical bytes.
"""

from __future__ import annotations

import os
import tempfile
import traceback

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from WatCon.scene import HIGHLY_CONSERVED, build_scene

#: Above this many structures the run is slow enough to warn about first.
SLOW_ABOVE = 60


def _chains_in_pdb(path):
    """Chain identifiers present in a PDB file, in order of first appearance."""
    seen = []
    try:
        with open(path, "r", errors="replace") as handle:
            for line in handle:
                if line.startswith("ENDMDL"):
                    break
                if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                    chain = line[21]
                    if chain not in seen and chain.strip():
                        seen.append(chain)
    except OSError:
        return []
    return seen


def _structure_names(directory):
    """Sorted stems of the structures in a directory, or [] if unreadable."""
    try:
        from WatCon.residue_index import list_structure_files

        names, _skipped = list_structure_files(directory)
    except Exception:                       # noqa: BLE001 - a bad path is normal
        return []
    return sorted(os.path.splitext(n)[0] for n in names)


class _Worker(QtCore.QThread):
    """Runs prepare + build_scene off the GUI thread."""

    progressed = QtCore.Signal(float, str)
    completed = QtCore.Signal(object)
    failed = QtCore.Signal(str)

    def __init__(self, options, parent=None):
        super().__init__(parent)
        self.options = options

    def run(self):                          # noqa: D102 - QThread entry point
        try:
            options = self.options
            structures = options["structures"]

            if options["prepare"]:
                from WatCon.prepare import prepare_directory

                self.progressed.emit(0.02, "Preparing structures...")
                prepared = os.path.join(options["out_dir"], "prepared")
                report = prepare_directory(
                    structures, prepared,
                    reference=options["reference"] or None,
                    chain_label=options["chain"],
                    min_identity=options["min_identity"],
                    verbose=False,
                )
                self.progressed.emit(
                    0.25, "Prepared %d of %d structures."
                    % (len(report.prepared), len(report.outcomes)))
                structures = prepared
                reference = None
            else:
                reference = options["reference"] or None

            def progress(fraction, message):
                # Preparation already used the first quarter of the bar.
                base = 0.25 if options["prepare"] else 0.0
                self.progressed.emit(base + fraction * (1.0 - base), message)

            scene = build_scene(
                structures, options["consurf"],
                out_dir=os.path.join(options["out_dir"], "view"),
                reference=reference,
                site_radius=options["site_radius"],
                min_cluster_samples=options["min_cluster_samples"],
                max_distance=options["max_distance"],
                highly_conserved=options["highly_conserved"],
                num_workers=options["num_workers"],
                progress=progress, verbose=False,
            )
            self.completed.emit(scene)
        except Exception:                   # noqa: BLE001 - report, never crash
            self.failed.emit(traceback.format_exc())


class WatConDialog(QtWidgets.QDialog):
    """Pick structures and parameters, run, and inspect the sites."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("WatCon + ConSurf")
        self.setMinimumWidth(760)
        self._worker = None
        self._scene = None
        self._build()

    # -- construction -------------------------------------------------------

    def _build(self):
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self._input_box())
        layout.addWidget(self._parameter_box())
        layout.addWidget(self._run_box())
        layout.addWidget(self._results_box(), 1)

    def _input_box(self):
        box = QtWidgets.QGroupBox("Input")
        form = QtWidgets.QFormLayout(box)

        self.structures_edit = QtWidgets.QLineEdit()
        self.structures_edit.setPlaceholderText(
            "folder of .pdb / .cif structures")
        self.structures_edit.editingFinished.connect(self._refresh_structures)
        browse = QtWidgets.QPushButton("Browse...")
        browse.clicked.connect(self._pick_structures)
        form.addRow("Structures", self._row(self.structures_edit, browse))

        self.ids_edit = QtWidgets.QLineEdit()
        self.ids_edit.setPlaceholderText("1AAX 7GSA   (downloads into the folder above)")
        fetch = QtWidgets.QPushButton("Fetch")
        fetch.clicked.connect(self._fetch)
        form.addRow("PDB ids", self._row(self.ids_edit, fetch))

        self.consurf_edit = QtWidgets.QLineEdit()
        self.consurf_edit.setPlaceholderText("*_consurf_grades.txt")
        consurf_browse = QtWidgets.QPushButton("Browse...")
        consurf_browse.clicked.connect(self._pick_consurf)
        form.addRow("ConSurf", self._row(self.consurf_edit, consurf_browse))

        self.reference_combo = QtWidgets.QComboBox()
        self.reference_combo.setEditable(True)
        self.reference_combo.currentTextChanged.connect(self._refresh_chains)
        self.chain_combo = QtWidgets.QComboBox()
        self.chain_combo.setEditable(True)
        self.chain_combo.addItem("A")
        pair = QtWidgets.QHBoxLayout()
        pair.setContentsMargins(0, 0, 0, 0)
        pair.addWidget(self.reference_combo, 3)
        pair.addWidget(QtWidgets.QLabel("chain"))
        pair.addWidget(self.chain_combo, 1)
        holder = QtWidgets.QWidget()
        holder.setLayout(pair)
        form.addRow("Reference", holder)

        self.prepare_check = QtWidgets.QCheckBox(
            "Prepare first: pick this chain, superpose, keep nearby waters")
        self.prepare_check.setChecked(True)
        self.prepare_check.setToolTip(
            "Uncheck only if the folder is already the output of "
            "`watcon prepare` -- every structure in one frame, one chain each.")
        form.addRow("", self.prepare_check)
        return box

    def _parameter_box(self):
        box = QtWidgets.QGroupBox("Parameters")
        grid = QtWidgets.QGridLayout(box)

        self.hbond_spin = self._spin(2.0, 5.0, 3.3, 0.1,
                                     "Donor-acceptor distance for a hydrogen "
                                     "bond, Angstrom.")
        self.radius_spin = self._spin(0.5, 5.0, 1.5, 0.1,
                                      "A water occupies a site within this "
                                      "distance of its centre, Angstrom.")
        self.min_samples_spin = QtWidgets.QSpinBox()
        self.min_samples_spin.setRange(2, 500)
        self.min_samples_spin.setValue(2)
        self.min_samples_spin.setToolTip(
            "Fewest waters that make a site. Raise it for large sets.")
        self.grade_spin = QtWidgets.QSpinBox()
        self.grade_spin.setRange(1, 9)
        self.grade_spin.setValue(HIGHLY_CONSERVED)
        self.grade_spin.setToolTip(
            "ConSurf grade at or above which a site counts as conserved. "
            "9 is most conserved.")
        self.workers_spin = QtWidgets.QSpinBox()
        self.workers_spin.setRange(1, 32)
        self.workers_spin.setValue(1)
        self.workers_spin.setToolTip("Parallel workers for network building.")
        self.identity_spin = self._spin(0.0, 1.0, 0.80, 0.05,
                                        "Reject a structure agreeing with the "
                                        "reference below this fraction.")

        for column, (label, widget) in enumerate([
            ("H-bond cutoff (A)", self.hbond_spin),
            ("Site radius (A)", self.radius_spin),
            ("Min waters / site", self.min_samples_spin),
        ]):
            grid.addWidget(QtWidgets.QLabel(label), 0, column * 2)
            grid.addWidget(widget, 0, column * 2 + 1)
        for column, (label, widget) in enumerate([
            ("Conserved grade >=", self.grade_spin),
            ("Min identity", self.identity_spin),
            ("Workers", self.workers_spin),
        ]):
            grid.addWidget(QtWidgets.QLabel(label), 1, column * 2)
            grid.addWidget(widget, 1, column * 2 + 1)
        return box

    def _run_box(self):
        box = QtWidgets.QWidget()
        row = QtWidgets.QHBoxLayout(box)
        row.setContentsMargins(0, 0, 0, 0)

        self.run_button = QtWidgets.QPushButton("Run")
        self.run_button.setDefault(True)
        self.run_button.clicked.connect(self._run)
        self.progress = QtWidgets.QProgressBar()
        self.progress.setRange(0, 100)
        self.status = QtWidgets.QLabel("Choose a folder of structures and a "
                                       "ConSurf file.")
        self.status.setWordWrap(True)

        row.addWidget(self.run_button)
        row.addWidget(self.progress, 1)
        outer = QtWidgets.QWidget()
        column = QtWidgets.QVBoxLayout(outer)
        column.setContentsMargins(0, 0, 0, 0)
        column.addWidget(box)
        column.addWidget(self.status)
        return outer

    def _results_box(self):
        box = QtWidgets.QGroupBox("Sites")
        column = QtWidgets.QVBoxLayout(box)

        self.summary = QtWidgets.QLabel("No results yet.")
        column.addWidget(self.summary)

        self.table = QtWidgets.QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(
            ["Site", "Grade", "Structures", "Waters", "Lining residues"])
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._focus_selected_site)
        column.addWidget(self.table, 1)

        self.conserved_only = QtWidgets.QCheckBox(
            "Only sites lined by a highly conserved residue")
        self.conserved_only.setChecked(True)
        self.conserved_only.stateChanged.connect(self._fill_table)
        column.addWidget(self.conserved_only)
        return box

    # -- small helpers ------------------------------------------------------

    @staticmethod
    def _row(*widgets):
        holder = QtWidgets.QWidget()
        row = QtWidgets.QHBoxLayout(holder)
        row.setContentsMargins(0, 0, 0, 0)
        for index, widget in enumerate(widgets):
            row.addWidget(widget, 1 if index == 0 else 0)
        return holder

    @staticmethod
    def _spin(low, high, value, step, tip):
        spin = QtWidgets.QDoubleSpinBox()
        spin.setRange(low, high)
        spin.setSingleStep(step)
        spin.setDecimals(2)
        spin.setValue(value)
        spin.setToolTip(tip)
        return spin

    def _say(self, message):
        self.status.setText(message)
        QtWidgets.QApplication.processEvents()

    # -- input ---------------------------------------------------------------

    def _pick_structures(self):
        chosen = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Folder of structures", self.structures_edit.text() or "")
        if chosen:
            self.structures_edit.setText(chosen)
            self._refresh_structures()

    def _pick_consurf(self):
        chosen, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "ConSurf grades file", self.consurf_edit.text() or "",
            "ConSurf grades (*grades*.txt);;Text files (*.txt);;All files (*)")
        if chosen:
            self.consurf_edit.setText(chosen)

    def _refresh_structures(self):
        """Fill the reference dropdown from whatever is in the folder."""
        directory = self.structures_edit.text().strip()
        names = _structure_names(directory) if directory else []
        current = self.reference_combo.currentText()
        self.reference_combo.blockSignals(True)
        self.reference_combo.clear()
        self.reference_combo.addItems(names)
        if current in names:
            self.reference_combo.setCurrentText(current)
        self.reference_combo.blockSignals(False)
        if names:
            self._say("%d structure(s) found. Reference: %s"
                      % (len(names), self.reference_combo.currentText()))
        elif directory:
            self._say("No structures in that folder. WatCon reads .pdb, .cif "
                      "and gzipped forms of either.")
        self._refresh_chains()

    def _refresh_chains(self):
        """Fill the chain dropdown from the reference structure itself."""
        directory = self.structures_edit.text().strip()
        reference = self.reference_combo.currentText().strip()
        if not directory or not reference:
            return
        path = None
        for suffix in (".pdb", ".ent", ".cif"):
            candidate = os.path.join(directory, reference + suffix)
            if os.path.isfile(candidate):
                path = candidate
                break
        chains = _chains_in_pdb(path) if path and path.endswith((".pdb", ".ent")) else []
        current = self.chain_combo.currentText()
        self.chain_combo.blockSignals(True)
        self.chain_combo.clear()
        self.chain_combo.addItems(chains or ["A"])
        if current in chains:
            self.chain_combo.setCurrentText(current)
        self.chain_combo.blockSignals(False)

    def _fetch(self):
        ids = self.ids_edit.text().replace(",", " ").split()
        if not ids:
            self._say("Type one or more PDB ids first, such as  1AAX 7GSA")
            return
        out_dir = self.structures_edit.text().strip()
        if not out_dir:
            out_dir = os.path.join(tempfile.gettempdir(), "watcon_structures")
            self.structures_edit.setText(out_dir)

        from WatCon.fetch import FetchError, fetch_structures

        self.run_button.setEnabled(False)
        try:
            self._say("Downloading %d structure(s)..." % len(ids))
            paths, failures = fetch_structures(ids, out_dir, verbose=False)
        except FetchError as error:
            self._say("Download failed. %s" % error)
            return
        finally:
            self.run_button.setEnabled(True)

        message = "Downloaded %d structure(s)." % len(paths)
        if failures:
            message += "  Could not fetch: %s" % ", ".join(i for i, _ in failures)
        self._say(message)
        self._refresh_structures()

    # -- running -------------------------------------------------------------

    def _options(self):
        structures = self.structures_edit.text().strip()
        consurf = self.consurf_edit.text().strip()
        if not structures or not os.path.isdir(structures):
            raise ValueError("Choose a folder of structures.")
        if not consurf or not os.path.isfile(consurf):
            raise ValueError("Choose a ConSurf *_consurf_grades.txt file. "
                             "WatCon does not contact the ConSurf server -- "
                             "run it yourself at consurf.tau.ac.il.")
        return {
            "structures": structures,
            "consurf": consurf,
            "reference": self.reference_combo.currentText().strip(),
            "chain": (self.chain_combo.currentText().strip() or "A")[:1],
            "prepare": self.prepare_check.isChecked(),
            "min_identity": self.identity_spin.value(),
            "max_distance": self.hbond_spin.value(),
            "site_radius": self.radius_spin.value(),
            "min_cluster_samples": self.min_samples_spin.value(),
            "highly_conserved": self.grade_spin.value(),
            "num_workers": self.workers_spin.value(),
            "out_dir": os.path.join(structures, "watcon_out"),
        }

    def _run(self):
        if self._worker is not None and self._worker.isRunning():
            self._say("Already running.")
            return
        try:
            options = self._options()
        except ValueError as error:
            self._say(str(error))
            return

        count = len(_structure_names(options["structures"]))
        if count > SLOW_ABOVE:
            answer = QtWidgets.QMessageBox.question(
                self, "This will take a while",
                "%d structures. Building the networks is the slow part and "
                "runs in the background -- PyMOL stays usable, but this can "
                "take several minutes.\n\nGo ahead?" % count,
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.Yes)
            if answer != QtWidgets.QMessageBox.Yes:
                return

        self.run_button.setEnabled(False)
        self.progress.setValue(0)
        self._say("Starting...")

        self._worker = _Worker(options, self)
        self._worker.progressed.connect(self._on_progress)
        self._worker.completed.connect(self._on_completed)
        self._worker.failed.connect(self._on_failed)
        self._worker.finished.connect(lambda: self.run_button.setEnabled(True))
        self._worker.start()

    def _on_progress(self, fraction, message):
        self.progress.setValue(int(max(0.0, min(1.0, fraction)) * 100))
        self._say(message)

    def _on_failed(self, detail):
        self.progress.setValue(0)
        last = [line for line in detail.strip().splitlines() if line.strip()]
        self._say("Failed: %s" % (last[-1] if last else "unknown error"))
        print(detail)                       # the full traceback, in PyMOL's log

    def _on_completed(self, scene):
        """Back on the GUI thread: display the scene and fill the table."""
        self._scene = scene
        self.progress.setValue(100)
        self._fill_table()
        try:
            cmd.delete("all")
            # Run the very file `watcon view` writes, so the plugin and the
            # command line cannot draw different pictures.
            cmd.do("@%s" % scene.pml_path.replace("\\", "/"))
        except Exception as error:          # noqa: BLE001
            self._say("Built the analysis, but PyMOL could not display it: %s"
                      % error)
            return
        self.summary.setText(scene.summary())
        self._say("Done. %s  Written to %s"
                  % (scene.summary(), scene.out_dir))

    # -- results -------------------------------------------------------------

    def _fill_table(self):
        scene = self._scene
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        if scene is None:
            self.table.setSortingEnabled(True)
            return

        sites = (scene.conserved_sites if self.conserved_only.isChecked()
                 else scene.sites)
        sites = sorted(sites, key=lambda s: (-s.grade, -s.occupancy_fraction))
        self.table.setRowCount(len(sites))
        for row, site in enumerate(sites):
            values = [
                site.cluster_id,
                site.grade if site.has_conservation else None,
                "%d/%d" % (site.n_structures_occupied, site.n_structures_total),
                site.occupancy,
                site.residue_label(),
            ]
            for column, value in enumerate(values):
                if value is None:
                    item = QtWidgets.QTableWidgetItem("no data")
                    item.setToolTip("No ConSurf score for any residue lining "
                                    "this site. That is not the same as low "
                                    "conservation.")
                elif isinstance(value, int):
                    item = QtWidgets.QTableWidgetItem()
                    item.setData(QtCore.Qt.DisplayRole, value)
                else:
                    item = QtWidgets.QTableWidgetItem(str(value))
                if column == 0:
                    item.setData(QtCore.Qt.UserRole, int(site.cluster_id))
                self.table.setItem(row, column, item)
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()
        self.summary.setText("%s  --  showing %d"
                             % (scene.summary(), len(sites)))

    def _focus_selected_site(self):
        """Fly to the selected site and show what lines it."""
        rows = self.table.selectionModel().selectedRows() if self.table.selectionModel() else []
        if not rows or self._scene is None:
            return
        item = self.table.item(rows[0].row(), 0)
        if item is None:
            return
        cluster_id = item.data(QtCore.Qt.UserRole)
        site = next((s for s in self._scene.sites
                     if int(s.cluster_id) == int(cluster_id)), None)
        if site is None:
            return

        # The table survives a `delete all`, so a row can be clicked when the
        # objects it refers to are gone. Say that, rather than letting PyMOL
        # report `Invalid selection name "sites"`, which reads like a bug.
        try:
            loaded = cmd.get_names("objects")
        except Exception:                   # noqa: BLE001
            loaded = []
        if "sites" not in loaded:
            self._say("Site %s is in the table, but the session is not loaded "
                      "-- press Run to draw it." % cluster_id)
            return

        selection = "sites and resi %d" % int(cluster_id)
        try:
            cmd.enable("sites")
            cmd.select("watcon_site", selection)
            parts = []
            for chain, resid, _icode, _grade in site.residues:
                parts.append("(chain %s and resi %d)" % (chain, resid)
                             if str(chain).strip() else "(resi %d)" % resid)
            if parts:
                cmd.select("watcon_lining",
                           "protein and (%s)" % " or ".join(parts))
                cmd.show("sticks", "watcon_lining and sidechain")
            cmd.zoom(selection, 8.0)
            cmd.deselect()
        except Exception as error:          # noqa: BLE001
            self._say("Could not focus site %s: %s" % (cluster_id, error))
            return

        grades = ", ".join(
            "%s%d%s" % (chain, resid, "" if grade is None else " (g%d)" % grade)
            for chain, resid, _icode, grade in site.residues)
        self._say("Site %s -- grade %s, occupied in %d/%d structures. Lined by: %s"
                  % (cluster_id,
                     site.grade if site.has_conservation else "no data",
                     site.n_structures_occupied, site.n_structures_total,
                     grades or "nothing recorded"))
