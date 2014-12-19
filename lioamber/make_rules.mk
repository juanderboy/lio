######################################################################
$(obj_path) :
	mkdir -p $@

$(obj_path)/%.o   :  %.f  |  $(obj_path)
	$(FC) -c $(FFLAGS) $< -o $@

$(obj_path)/%.mod :  %.f  |  $(obj_path)
	$(FC) -c $(FFLAGS) $< -o $(@:%.mod=%.o)

.PHONY: clean
clean:
	rm -rf liblio-g2g.so $(obj_path)
######################################################################
