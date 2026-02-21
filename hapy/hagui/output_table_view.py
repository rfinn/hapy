
from astroyp.table import Table

class OutputTableView():
    """ methods that interact with the gui  """
    
    def update_gui_table(self):

        self.ui.tableWidget.setColumnCount(len(self.table.columns))
        self.ui.tableWidget.setRowCount(len(self.table))
        
        self.ui.tableWidget.setHorizontalHeaderLabels(self.table.colnames)

        if self.igal is not None:
            #print(self.ui.commentLineEdit.text())
            self.table['COMMENT'][self.igal] = str(self.ui.commentLineEdit.text())
            
        for col, c in enumerate(self.table.columns):
            #item = self.ui.tableWidget.horizontalHeaderItem(col)
            #item.setText(_translate("MainWindow", self.table.columns[col].name))
            for row in range(len(self.table[c])):
                item = self.table[row][col]
                self.ui.tableWidget.setItem(row,col,QtWidgets.QTableWidgetItem(str(item)))

    def update_gui_table_cell(self,row,col,item):
        # row will be igal, so that's easy
        # how to make it easy to get the right column number?
        #
        # easiest is to pass the column name, and
        # then match the column name to find the

        colmatch = False
        for i,c in enumerate(self.table.colnames):
            if c == col:
                ncol = i
                colmatch = True
                break
        if colmatch:    
            self.ui.tableWidget.setItem(row,ncol,QtWidgets.QTableWidgetItem(str(item)))
            self.table[row][col]=item
        else:
            print('could not match column name ',col)
        # right now, write_fits_table calls update_gui_table_cell!!!
        # this this is recursive :(
        #self.write_fits_table()
