ALL: quad_x310_logger.py quad_logger_display.py quad_x310_logger_fft_logger.py

quad_x310_logger.py: quad_x310_logger.grc
	grcc quad_x310_logger.grc

quad_x310_logger_fft_logger.py: quad_x310_logger.grc
	grcc quad_x310_logger.grc

quad_logger_display.py: quad_logger_display.grc
	grcc quad_logger_display.grc

install: quad_x310_logger.py quad_logger_display.py quad_x310_logger_fft_logger.py
	cp quad_x310_logger.py /usr/local/bin
	chmod 755 /usr/local/bin/quad_x310_logger.py
	cp quad_x310_logger_fft_logger.py /usr/local/bin/
	cp quad_logger_display.py /usr/local/bin
	chmod 755 /usr/local/bin/quad_logger_display.py

clean:
	rm -f quad_x310_logger.py
	rm -f quad_logger_display.py
	rm -f quad_x310_logger_fft_logger.py
