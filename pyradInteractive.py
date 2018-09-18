import pyrad


class Menu:
    def __index__(self, title, *entries):
        self.title = title
        self.entries = entries

    def displayMenu(self):
        print('***\tPyrad v%s\t\t\t***', pyrad.VERSION)
        i = 1
        for entry in self.entries:
            print('%s.  %s' % (i, entry.text))
        print('0.  Exit')

    def receiveInput(self):

