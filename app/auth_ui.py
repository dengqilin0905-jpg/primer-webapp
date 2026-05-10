import time, os, base64, hashlib
import streamlit as st
from db import get_conn, get_app_setting, set_app_setting

def _hash_password(pw: str, salt: bytes) -> str:
    dk = hashlib.pbkdf2_hmac("sha256", pw.encode("utf-8"), salt, 120_000)
    return base64.b64encode(salt + dk).decode("utf-8")

def _verify_password(pw: str, stored: str) -> bool:
    raw = base64.b64decode(stored.encode("utf-8"))
    salt, dk = raw[:16], raw[16:]
    dk2 = hashlib.pbkdf2_hmac("sha256", pw.encode("utf-8"), salt, 120_000)
    return dk2 == dk

def get_current_user():
    if "user" in st.session_state:
        return st.session_state["user"]
    return None

def logout():
    st.session_state.pop("user", None)
    st.session_state.pop("mapping_cart", None)

def login_box():
    st.sidebar.header("Account")
    tab1, tab2 = st.sidebar.tabs(["Login", "Register"])
    if "remember_username" not in st.session_state:
        conn = get_conn()
        remembered = get_app_setting(conn, "remember_username", "")
        conn.close()
        st.session_state["remember_username"] = remembered or ""

    with tab1:
        if "login_u" not in st.session_state and st.session_state["remember_username"]:
            st.session_state["login_u"] = st.session_state["remember_username"]
        u = st.text_input("Username", key="login_u")
        p = st.text_input("Password", type="password", key="login_p")
        remember_u = st.checkbox(
            "Remember username on this device",
            value=bool(st.session_state.get("remember_username")),
            key="remember_username_check",
        )
        if st.button("Login"):
            conn = get_conn()
            row = conn.execute("SELECT * FROM users WHERE username=?", (u.strip(),)).fetchone()
            if not row:
                st.error("User not found.")
            else:
                if _verify_password(p, row["password_hash"]):
                    st.session_state["user"] = {"id": row["id"], "username": row["username"]}
                    if remember_u and u.strip():
                        set_app_setting(conn, "remember_username", u.strip())
                        st.session_state["remember_username"] = u.strip()
                    else:
                        set_app_setting(conn, "remember_username", None)
                        st.session_state["remember_username"] = ""
                    st.success("Logged in.")
                    st.rerun()
                else:
                    st.error("Wrong password.")
            conn.close()

    with tab2:
        u = st.text_input("New username", key="reg_u")
        p1 = st.text_input("New password", type="password", key="reg_p1")
        p2 = st.text_input("Repeat password", type="password", key="reg_p2")
        if st.button("Create account"):
            if not u.strip():
                st.error("Username required.")
            elif p1 != p2 or not p1:
                st.error("Passwords do not match.")
            else:
                salt = os.urandom(16)
                ph = _hash_password(p1, salt)
                conn = get_conn()
                try:
                    conn.execute(
                        "INSERT INTO users(username, password_hash, created_at) VALUES(?,?,?)",
                        (u.strip(), ph, time.strftime("%Y-%m-%d %H:%M:%S")),
                    )
                    conn.commit()
                    st.success("Account created. Please login.")
                except Exception as e:
                    st.error(f"Create failed: {e}")
                finally:
                    conn.close()
