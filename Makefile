include gmsl

RUNDIR := runs
INPUTDIR := $(RUNDIR)/input
DATADIR := $(RUNDIR)/data
OUTPUTDIR := $(RUNDIR)/output
REPORTDIR := $(RUNDIR)/reports-pub
RESULTSDIR := $(RUNDIR)/results
CODEDIR := code

HUMAN_PREFIX := 9606
YEAST_PREFIX := 559292
PPITRIM_SUFFIX := mitab.2012-11-23.txt.ppiTrim.txt
BETAS := 0.0 4.0 5.0 2.5 6.5
ENRICHED_BETAS := 4.0

SADDLESUM := /home/stojmira/virtualenvs/qmbp20120920/bin

gruns_out = $(shell python -c "print '{}_{:03d}_{:03d}'.format('$(1)',int($(2)*100),int($(3)*100))")
supp_out = $(shell python -c "print '{}_supp_{:03d}'.format('$(1)',int($(2)*100))")

havugimana := $(INPUTDIR)/havugimana2012.ppiTrim.txt
human_db := $(INPUTDIR)/$(HUMAN_PREFIX).$(PPITRIM_SUFFIX)
yeast_db := $(INPUTDIR)/$(YEAST_PREFIX).$(PPITRIM_SUFFIX)
psi_mi_obo := $(INPUTDIR)/psi-mi.obo
yeast_etd := $(INPUTDIR)/Saccharomyces_cerevisiae.etd
human_etd := $(INPUTDIR)/Homo_sapiens.etd
human_db_full := $(DATADIR)/$(HUMAN_PREFIX).fulldb.txt
yeast_db_full := $(DATADIR)/$(YEAST_PREFIX).fulldb.txt

species_prefixes := $(HUMAN_PREFIX) $(YEAST_PREFIX)
yeast_runs := $(foreach alpha,0.0 0.15,$(call supp_out,yeast,$(alpha))) $(foreach beta,$(BETAS),$(call gruns_out,yeast,$(beta),0.15))
human_runs := $(foreach alpha,0.0 0.15,$(call supp_out,human,$(alpha))) $(foreach beta,$(BETAS),$(call gruns_out,human,$(beta),0.15))

summary_dags := supp_000 supp_015 000_015 400_015
summary_lbls := Z P A N4
yeast_summary_dags := $(foreach sfx,$(summary_dags), $(OUTPUTDIR)/yeast_$(sfx).dagf)
human_summary_dags := $(foreach sfx,$(summary_dags), $(OUTPUTDIR)/human_$(sfx).dagf)
yeast_summary_lbls := $(foreach sfx,$(summary_lbls), sce$(sfx))
human_summary_lbls := $(foreach sfx,$(summary_lbls), hsa$(sfx))
summary_param = '$(1)@$(2)'
yeast_enriched := $(foreach beta,$(ENRICHED_BETAS),$(call gruns_out,yeast,$(beta),0.15))
human_enriched := $(foreach beta,$(ENRICHED_BETAS),$(call gruns_out,human,$(beta),0.15))


yeast_dags := $(foreach item,$(yeast_runs),$(OUTPUTDIR)/$(item).dagf)
human_dags := $(foreach item,$(human_runs),$(OUTPUTDIR)/$(item).dagf)

all_cmembers := $(foreach spc,$(species_prefixes),$(DATADIR)/$(spc).cmembers.pkl)
all_pplinks := $(foreach spc,$(species_prefixes),$(DATADIR)/$(spc).pplinks.txt)


# ------------------------------------------------------------------------------------
# Templates for custom recipies
# ------------------------------------------------------------------------------------

define GRUNS_template
# $(1) - species tag
# $(2) - beta
# $(3) - alpha
# $(4) - species taxid
$(OUTPUTDIR)/$(call gruns_out,$(1),$(2),$(3)).results.pkl: $(DATADIR)/$(4).cmembers.pkl $(DATADIR)/$(4).pplinks.txt | $(OUTPUTDIR)
	$(CODEDIR)/Gruns.py -p 0 -c 0 $$@ $$+ $(3) $(2) 0.0 > $$@.log
endef

define SUPPORT_template
# $(1) - species tag
# $(2) - alpha
# $(3) - species taxid
$(OUTPUTDIR)/$(call supp_out,$(1),$(2)).results.pkl: $(DATADIR)/$(3).cmembers.pkl | $(OUTPUTDIR)
	$(CODEDIR)/plainmerge.py $$@ $$< $(2) 0.0 > $$@.log
endef

define ENRICH_DAG_template
# $(1) - run name
# $(2) - term database
$(OUTPUTDIR)/$(1).enrich: $(OUTPUTDIR)/$(1).dagf
	$(CODEDIR)/dag_enrich.py -v -x $(SADDLESUM) $$< $(2)
endef


.PHONY: all datasets cmembers pplinks yeast_dags human_dags
.PHONY: reports adb_reports yeast_reports human_reports
.PHONY: clean


all: reports

# ------------------------------------------------------------------------------------
# Recipies for generating datasets
# ------------------------------------------------------------------------------------

datasets: $(human_db_full) $(yeast_db_full)

cmembers: $(all_cmembers)

pplinks: $(all_pplinks)



$(human_db_full): $(human_db) $(havugimana) | $(DATADIR)
	tail -n +2 $(havugimana) | cat $(human_db) - > $(human_db_full)

$(yeast_db_full): $(yeast_db) | $(DATADIR)
	cp $(realpath $(yeast_db)) $(yeast_db_full)

$(DATADIR)/%.cmembers.pkl: $(DATADIR)/%.fulldb.txt
	$(CODEDIR)/adb_create.py $< $@

$(DATADIR)/%.pplinks.txt: $(DATADIR)/%.fulldb.txt
	$(CODEDIR)/pplinks.py $< $(psi_mi_obo) $@


# ------------------------------------------------------------------------------------
# Recipies for generating DAGs (iterative merge, filtering, enrichment)
# ------------------------------------------------------------------------------------

yeast_dags: $(yeast_dags) $(foreach item,$(yeast_enriched), $(OUTPUTDIR)/$(item).enrich)

human_dags: $(human_dags) $(foreach item,$(human_enriched), $(OUTPUTDIR)/$(item).enrich)


$(foreach beta,$(BETAS),$(eval $(call GRUNS_template,yeast,$(beta),0.15,$(YEAST_PREFIX))))
$(foreach beta,$(BETAS),$(eval $(call GRUNS_template,human,$(beta),0.15,$(HUMAN_PREFIX))))
$(foreach alpha,0.0 0.15,$(eval $(call SUPPORT_template,yeast,$(alpha),$(YEAST_PREFIX))))
$(foreach alpha,0.0 0.15,$(eval $(call SUPPORT_template,human,$(alpha),$(HUMAN_PREFIX))))


$(OUTPUTDIR)/%.dagf: $(OUTPUTDIR)/%.results.pkl
	$(CODEDIR)/dag_filter.py -p 2.5 $< $@

$(foreach item,$(yeast_enriched),$(eval $(call ENRICH_DAG_template,$(item),$(yeast_etd))))
$(foreach item,$(human_enriched),$(eval $(call ENRICH_DAG_template,$(item),$(human_etd))))


# ------------------------------------------------------------------------------------
# Recipies for generating reports
# ------------------------------------------------------------------------------------

reports: adb_reports yeast_reports human_reports

adb_reports: $(foreach spc,$(species_prefixes),$(REPORTDIR)/$(spc).p2c.txt) $(foreach spc,$(species_prefixes),$(REPORTDIR)/$(spc).p2p.txt) $(foreach spc,$(species_prefixes),$(REPORTDIR)/$(spc).raw_complexes.txt)

yeast_reports: $(foreach item,$(yeast_runs),$(REPORTDIR)/$(item).node_details.txt $(REPORTDIR)/$(item).hier.final.edges.tab $(REPORTDIR)/$(item).hier.expanded.edges.tab) $(foreach item,$(yeast_enriched), $(REPORTDIR)/$(item).node_enrich.avgw.txt $(REPORTDIR)/$(item).node_enrich.wsum.txt $(REPORTDIR)/$(item).component_enrich.avgw.txt $(REPORTDIR)/$(item).component_details.txt)

human_reports: $(foreach item,$(human_runs),$(REPORTDIR)/$(item).node_details.txt $(REPORTDIR)/$(item).hier.final.edges.tab $(REPORTDIR)/$(item).hier.expanded.edges.tab) $(foreach item,$(human_enriched), $(REPORTDIR)/$(item).node_enrich.avgw.txt $(REPORTDIR)/$(item).node_enrich.wsum.txt $(REPORTDIR)/$(item).component_enrich.avgw.txt $(REPORTDIR)/$(item).component_details.txt)



$(REPORTDIR)/%.raw_complexes.txt: $(DATADIR)/%.cmembers.pkl | $(REPORTDIR)
	$(CODEDIR)/adb_report.py details $< $@

$(REPORTDIR)/%.p2c.txt: $(DATADIR)/%.cmembers.pkl | $(REPORTDIR)
	$(CODEDIR)/adb_report.py p2c $< $@

$(REPORTDIR)/%.p2p.txt: $(DATADIR)/%.pplinks.txt | $(REPORTDIR)
	$(CODEDIR)/adb_report.py p2p $< $@

$(REPORTDIR)/%.node_details.txt: $(OUTPUTDIR)/%.dagf | $(REPORTDIR)
	$(CODEDIR)/dag_report.py node_details $< $@

$(REPORTDIR)/%.component_details.txt: $(OUTPUTDIR)/%.dagf | $(REPORTDIR)
	$(CODEDIR)/dag_report.py original_component_details $< $@

$(REPORTDIR)/%.node_enrich.avgw.txt: $(OUTPUTDIR)/%.dagf $(OUTPUTDIR)/%.enrich | $(REPORTDIR)
	$(CODEDIR)/dag_report.py node_enrichment $< $@ avgw

$(REPORTDIR)/%.node_enrich.wsum.txt: $(OUTPUTDIR)/%.dagf $(OUTPUTDIR)/%.enrich | $(REPORTDIR)
	$(CODEDIR)/dag_report.py node_enrichment $< $@ wsum

$(REPORTDIR)/%.component_enrich.avgw.txt: $(OUTPUTDIR)/%.dagf $(OUTPUTDIR)/%.enrich | $(REPORTDIR)
	$(CODEDIR)/dag_report.py original_component_enrichment $< $@ avgw

$(REPORTDIR)/%.hier.final.edges.tab: $(OUTPUTDIR)/%.dagf | $(REPORTDIR)
	$(CODEDIR)/dag_report.py final_dag_data $< - $(subst .edges.tab,,$@)

$(REPORTDIR)/%.hier.expanded.edges.tab: $(OUTPUTDIR)/%.dagf | $(REPORTDIR)
	$(CODEDIR)/dag_report.py expanded_dag_data $< $@



# ------------------------------------------------------------------------------------
# Other recipies
# ------------------------------------------------------------------------------------

$(DATADIR):
	mkdir -p $(DATADIR)

$(OUTPUTDIR):
	mkdir -p $(OUTPUTDIR)

$(REPORTDIR):
	mkdir -p $(REPORTDIR)

$(RESULTSDIR):
	mkdir -p $(RESULTSDIR)

clean:
	-rm -rf $(DATADIR)
	-rm -rf $(OUTPUTDIR)
	-rm -rf $(REPORTDIR)
	-rm -rf $(RESULTSDIR)
