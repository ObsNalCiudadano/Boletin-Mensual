# Boletin-Mensual
Boletín mensual del Observatorio Nacional Ciudadano

Aquí se encuentran el código que el Observatorio Nacional Ciudadano utiliza para generar su boletín mensual de seguridad en México.

El archivo boletin_setup.R se corre primero para generar las visualizaciones y subirlas al servidor ya que el correo electrónico utiliza enlaces desde el html.

Los archivos de las tablas no se corren individualmente pero el primer archivo los usa.

Al finalizar el primer script se puede correr render_boletin.R para generar el código html que se usa en el correo del boletín mensual.
