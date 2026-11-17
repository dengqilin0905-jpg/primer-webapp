"""Generate empty generic software templates; no provider customer records."""
from pathlib import Path
from openpyxl import Workbook
from openpyxl.workbook.defined_name import DefinedName
from openpyxl.worksheet.datavalidation import DataValidation
from openpyxl.styles import Font,PatternFill

def generate(destination):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    wb=Workbook();ws=wb.active;ws.title='寡核苷酸合成订购单'
    ws['A1']='Generic primer order template; confirm current provider requirements before submission.'
    for row,label in enumerate(['Orderer','Mobile','Email','Legal person','Phone','Institution','Invoice','Payment','Invoice delivery','Shipping address','Note'],4):ws.cell(row,1,label);ws.cell(row,2,0)
    headers=['ID','Primer名称','序列','碱基数','分装管数','提供总量','计量单位','纯化方式',"5'端修饰","3'端修饰",'双标记修饰','类型','备注']
    for col,value in enumerate(headers,1):ws.cell(17,col,value)
    for row in range(18,224):ws.cell(row,4,f'=LEN(C{row})')
    for col in ['A','B','C','D','E','F','G','H','I','J','K','L','M']:ws.column_dimensions[col].width=20
    wb.properties.creator='Primer WebApp';wb.properties.lastModifiedBy='Primer WebApp';wb.save(destination/'order_generic.xlsx')
    wb=Workbook();ws=wb.active;ws.title='上传模版';hidden=wb.create_sheet('请勿删除！');hidden.sheet_state='hidden'
    ws['A1']='Generic sequencing template; configure the current company workbook for formal orders.'
    headers=['* 样品名称','* 样品类型','* 测序项目','* 参考序列或引物序列','* 质粒总长度','* 片段长度','* 抗性','是否返样','返样类型','备注']
    for col,value in enumerate(headers,1):ws.cell(3,col,value)
    for row,value in enumerate(['氨苄','卡那','氯霉素','壮观霉素'],1):hidden.cell(row,1,value)
    wb.defined_names.add(DefinedName('Antibiotic',attr_text="'请勿删除！'!$A$1:$A$4"))
    for col,formula in [('B','"平板,菌液,质粒"'),('C','"插入片段测通,全质粒测序"'),('G','Antibiotic'),('H','"是,否"'),('I','"质粒+菌液,质粒,菌液"')]:
        dv=DataValidation(type='list',formula1=formula);dv.errorTitle='Invalid option';dv.error='Select a listed value';dv.showErrorMessage=True;ws.add_data_validation(dv);dv.add(col+'4:'+col+'60')
    for col in 'ABCDEFGHIJ':ws.column_dimensions[col].width=24
    wb.properties.creator='Primer WebApp';wb.properties.lastModifiedBy='Primer WebApp';wb.save(destination/'sequencing_generic.xlsx')
    return destination
if __name__=='__main__':generate(Path(__file__).resolve().parents[1]/'resources/templates')
