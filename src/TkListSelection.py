#!/usr/bin/python

from Tkinter import *
import os

class TkListSelectionClass:
    def __init__(self, items, title):
        self.items = items
        self.selection = None
        master = Tk()

        title = Label(master, text=title)

        frame = Frame(master)
        scrollbar = Scrollbar(frame, orient=VERTICAL)
        listbox = Listbox(frame, selectmode=SINGLE, yscrollcommand=scrollbar.set)
        scrollbar.config(command=listbox.yview)
        scrollbar.pack(side=RIGHT, fill=Y)

        for item in items:
            listbox.insert(END, item)

        listbox.selection_set(0)
        listbox.bind("<Double-Button-1>", self.callback_ok)
        listbox.bind("<Return>", self.callback_ok)
        listbox.bind("<Cancel>", self.callback_cancel)
        ok_button = Button(master, text="OK", command=self.ok)

        title.pack()
        frame.pack()
        listbox.pack()
        ok_button.pack()

        self.master = master
        self.listbox = listbox
        mainloop()

    def get_selection(self):
        index = int(self.listbox.curselection()[0])
        return self.items[index]

    def ok(self):
        self.selection = self.get_selection()
        self.master.destroy()
        
    def callback_ok(self, event):
        self.ok()
        
    def callback_cancel(self, event):
        self.selection = None
        self.master.destroy()

def TkListSelection(items, title="Choose an item from the list:"):
    list_selection = TkListSelectionClass(items, title)
    return list_selection.selection
    
def test():
    while (True):
        choice = TkListSelection(["a","b","c"])
        if (choice != None):
            print choice
        else:
            break

if __name__ == '__main__': test()
