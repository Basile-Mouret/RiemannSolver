TYPST = typst 
FILENAME = Notes_on_Toro

$(FILENAME).pdf : $(FILENAME).typ chapters/*
	$(TYPST) compile $(FILENAME).typ $@

clean: 
	@rm -f $(FILENAME).pdf 
