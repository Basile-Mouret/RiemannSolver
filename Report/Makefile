TYPST=typst 
FILENAME = report_internship_cagire

$(FILENAME).pdf : $(FILENAME).typ chapters/*
	$(TYPST) compile $(FILENAME).typ $@

clean: 
	@rm -f $(FILENAME).pdf 
