from rich import print as rprint
from rich.panel import Panel
from rich.align import Align

class Messages:
    """ General stuff to print """

    def hello():
        w_message = "Ephemeris Fix"
        w_panel = Panel(w_message, padding=2)
        w_panel = Align.left(w_panel, vertical="middle")
        rprint(w_panel)

    def goodbye():
        bye_message = "Session Complete"
        bye_panel = Panel(bye_message, padding=2)
        bye_panel = Align.left(bye_panel, vertical="middle")
        rprint(bye_panel)