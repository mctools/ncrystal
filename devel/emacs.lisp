;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Loading this file from an active emacs session with M-x load-file will
;; provide a few hopefully useful developer utilities.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Make sure Meta+O key can be used to jump between files in ncrystal repo
;; (jumps to corresponding .cc/.hh file unless cursor is on a line with an
;; include statement):

(setq ncrystal-path-emacslispfile (or load-file-name (or buffer-file-name "SELFFILENOTFOUND")))
(setq ncrystal-path-develdir (file-name-parent-directory ncrystal-path-emacslispfile ))
(setq ncrystal-path-reporoot (file-name-parent-directory ncrystal-path-develdir ))
(setq ncrystal-path-coredir (format "%s/ncrystal_core" ncrystal-path-reporoot))
(setq cc-search-directories (list (format "%s/include" ncrystal-path-coredir)
                                  (format "%s/include/NCrystal/internal/*" ncrystal-path-coredir)
                                  (format "%s/include/NCrystal/*" ncrystal-path-coredir)
                                  "."
                                  (format "%s/src/*" ncrystal-path-coredir)))
(global-set-key "\M-o" 'ff-find-other-file)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Major mode for .ncmat files:
(setq ncmat-sections '("CELL" "ATOMPOSITIONS" "SPACEGROUP" "DEBYETEMPERATURE" "DYNINFO" "TEMPERATURE"
                       "DENSITY" "ATOMDB" "STATEOFMATTER" "OTHERPHASES" "CUSTOM_[A-Z]+"))
(setq ncmat-fields '("lengths" "angles" "cubic" "element" "fraction" "type" "debye_temp" "temperature"
                     "sab_scaled" "sab" "alphagrid" "betagrid" "egrid" "vdos_egrid" "vdos_density"))
(setq ncmat-highlights
      `(("^NCMAT\s*v[0-9]+[a-z0-9]*" . font-lock-function-name-face)
        (,(concat "^@\\(" (mapconcat 'identity ncmat-sections "\\|") "\\)") . font-lock-keyword-face)
        (,(concat "^\s*\\(" (mapconcat 'identity ncmat-fields "\\|") "\\)") . font-lock-constant-face)
        ))
(setq ncmat-mode-syntax-table
      (let ( (synTable (make-syntax-table)))
        (modify-syntax-entry ?# "<" synTable)
        (modify-syntax-entry ?\n ">" synTable)
        synTable)
      )
(define-derived-mode ncmat-mode fundamental-mode "NCMAT"
  "major mode for editing .ncmat files."
  (setq font-lock-defaults '(ncmat-highlights));keywords
  (set-syntax-table ncmat-mode-syntax-table);comments
  )
(add-to-list 'auto-mode-alist '("\\.ncmat$" . ncmat-mode))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
