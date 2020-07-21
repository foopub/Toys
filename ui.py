"""
Chat UI project for learning the curses library.
"""

from _curses import window
import curses
from curses.textpad import Textbox
import asyncio

win = {}

def colours() -> None:
    """
    Init all the colours and colour pairs. Can only be called after
    initscr(), hence inside the wrapped main!
    """
    #terminal default               #white on black (0)
    curses.init_pair(1, 1, 8)       #red on brblack (1)
    curses.init_pair(2, 127, 233)   #violet on gray (2)
    curses.init_pair(3, 9, 234)     #orange on gray (3)

    scolour = curses.color_pair(3)
    tcolour = curses.color_pair(1)
    mcolour = curses.color_pair(2)

    win['s'].attrset(scolour)     #Set the colour for text
    win['s'].bkgd(' ', scolour)   #Set the background 
    win['m'].attrset(mcolour)     #Set the colour for text
    win['m'].bkgd(' ', mcolour)   #Set the background 
    win['t'].attrset(tcolour)     #Set the colour for text
    win['t'].bkgd(' ', tcolour)   #Set the background 

def draw_boxes(stdscr: window):
    """
    This method draws the boxes on startup
    """

    h = curses.LINES    #Total height 
    w = curses.COLS     #Total width

    #Some debug values, remove these later
    text = "This is the start of a new window"

    # Clear screen
    stdscr.clear()
    stdscr.addstr(0, 0, 'This is the start of a new window')
    stdscr.noutrefresh()

    theight = 2
    swidth = 15
    mheight = h-theight-2
    mwidth = w-swidth-2

    win['t'] = curses.newwin(theight, mwidth, h-1-theight, 1)
    win['t'].addstr(0, 0, "This is the messages windows")
    
    win['m'] = curses.newwin(mheight, mwidth, 1, 1)
    win['m'].addstr(0, 0, "Type something here.")
    win['m'].timeout(0)

    win['s'] = curses.newwin(h-2, swidth, 1,  w-swidth-1)
    win['s'].addstr(0, 0, "This is a side pane")

async def textbox(handler):
    """
    Handles user input to the textbox.
    """
    while True:
        win['t'].erase()
        box = Textbox(win['t'])
        box.edit()
        message = box.gather()
        await handler(message)

async def messages(input):
    """
    Handles input to the messages box.
    """
    while True:
        message = await input()
        win['m'].addstr(message)
        win['m'].noutrefresh()

#Example handlers for testing----------------------------------
from socket import socketpair
rsock, wsock = socketpair()

async def writeToSock(message: str):
    """
    Example output handler.
    """
    wsock.send(message.encode())
    await asyncio.sleep(0)

async def readFromSock() -> str:
    """
    Example imput handler.
    """
    await asyncio.sleep(0)
    data = rsock.recv(100)
    return data.decode()
#--------------------------------------------------------------

async def main() ->  None:
    """
    Example main loop for testing.
    """
    loop = asyncio.get_running_loop()

    await asyncio.gather(messages(readFromSock) , textbox(writeToSock))

def init(main):
    #This is based on the curses.wrapper() func
    try:
    #Initialise curses
        stdscr = curses.initscr()
        curses.noecho()
        curses.cbreak()
        stdscr.keypad(True)
        #Do things here
        draw_boxes(stdscr)
        try:
            curses.start_color()
            colours()
        except:
            pass
        for i in win.values():
            i.noutrefresh()
        curses.doupdate()

        asyncio.run(main())
        
    finally:
    #Terminate curses
        if "stdscr" in locals():
            stdscr.keypad(0)
            curses.echo()
            curses.nocbreak()
            curses.endwin()

if __name__ == "__main__":
    init(main)
