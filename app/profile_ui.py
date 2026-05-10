import time
import streamlit as st
from db import get_conn

PAYMENT_OPTIONS = ["银行转账", "协议支付", "专项经费", "货到付款"]
YESNO_OPTIONS = ["是", "否"]

def profile_page(user):
    st.subheader("Profile (Template header auto-fill)")
    conn = get_conn()
    row = conn.execute("SELECT * FROM user_profiles WHERE user_id=?", (user["id"],)).fetchone()
    prof = dict(row) if row else {}

    col1, col2 = st.columns(2)
    with col1:
        orderer_name = st.text_input("*订购人姓名 (B4)", value=prof.get("orderer_name") or "")
        mobile = st.text_input("*订购人手机 (B5)", value=prof.get("mobile") or "")
        email = st.text_input("*订购人E-MAIL (B6)", value=prof.get("email") or "")
        legal_person = st.text_input("*负责人/法人 (B7)", value=prof.get("legal_person") or "")
        phone = st.text_input("*固定电话 (B8)", value=prof.get("phone") or "")
    with col2:
        institution = st.text_input("*订购人单位 (B9)", value=prof.get("institution") or "")
        invoice_title = st.text_input("*发票抬头 (B10)", value=prof.get("invoice_title") or "")
        payment_method = st.selectbox("*付款方式 (B11)", PAYMENT_OPTIONS,
                                      index=PAYMENT_OPTIONS.index(prof.get("payment_method")) if prof.get("payment_method") in PAYMENT_OPTIONS else 0)
        invoice_with_goods = st.selectbox("随货开票 (B12)", YESNO_OPTIONS,
                                          index=YESNO_OPTIONS.index(prof.get("invoice_with_goods")) if prof.get("invoice_with_goods") in YESNO_OPTIONS else 1)
        shipping_address = st.text_area("*收货地址 (B13)", value=prof.get("shipping_address") or "", height=90)
        header_note = st.text_area("表头备注 (B14)", value=prof.get("header_note") or "", height=70)

    if st.button("Save profile"):
        now = time.strftime("%Y-%m-%d %H:%M:%S")
        conn.execute(
            '''
            INSERT INTO user_profiles(user_id, orderer_name, mobile, email, legal_person, phone,
                                      institution, invoice_title, payment_method, invoice_with_goods,
                                      shipping_address, header_note, updated_at)
            VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)
            ON CONFLICT(user_id) DO UPDATE SET
                orderer_name=excluded.orderer_name,
                mobile=excluded.mobile,
                email=excluded.email,
                legal_person=excluded.legal_person,
                phone=excluded.phone,
                institution=excluded.institution,
                invoice_title=excluded.invoice_title,
                payment_method=excluded.payment_method,
                invoice_with_goods=excluded.invoice_with_goods,
                shipping_address=excluded.shipping_address,
                header_note=excluded.header_note,
                updated_at=excluded.updated_at
            ''',
            (
                user["id"],
                orderer_name.strip() or None,
                mobile.strip() or None,
                email.strip() or None,
                legal_person.strip() or None,
                phone.strip() or None,
                institution.strip() or None,
                invoice_title.strip() or None,
                payment_method,
                invoice_with_goods,
                shipping_address.strip() or None,
                header_note.strip() or None,
                now,
            ),
        )
        conn.commit()
        st.success("Saved.")
    conn.close()
