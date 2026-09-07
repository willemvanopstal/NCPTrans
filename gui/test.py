from nicegui import ui
ui.label('Hellow')
ui.run(host='0.0.0.0', port=8080, forwarded_allow_ips='*')
